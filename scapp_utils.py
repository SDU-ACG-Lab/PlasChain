from __future__ import division
import numpy as np
import math,copy,os
import networkx as nx
import re, pysam,sys
import logging
from plasclass import plasclass
import multiprocessing as mp
from multiprocessing import Manager
from heapq import heappush, heappop
from itertools import count
import time
from concurrent.futures import ThreadPoolExecutor
from parse_plasmid_scores import transformByLength
from sklearn.metrics.pairwise import cosine_similarity
from collections import defaultdict, Counter, deque
from node_feature import *
from functools import wraps

import PARAMS

complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
logger = logging.getLogger("scapp_logger")


worker_data = {}

def init_worker(path_dict, node_score_dict, node_vec_dict, node_support_dict, SEQS, G):
    """
    这个函数会在每个 Worker 启动时执行一次。
    它接收主进程传来的最新 G，并存入自己的全局变量。
    """
    global worker_data
    worker_data['G'] = G
    worker_data['SEQS'] = SEQS
    worker_data['scores'] = node_score_dict
    worker_data['path'] = path_dict
    worker_data['vec'] = node_vec_dict
    worker_data['support'] = node_support_dict

def readfq(fp): # this is a generator function
    """ # lh3's fast fastX reader:
        https://github.com/lh3/readfq/blob/master/readfq.py
    """
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def get_node_freq_vec(G,SEQs):
    """ Annotate each node in the graph with its frequency vector
    """
    node_vec_dict = {}
    unique_canonical_nodes = {canonicalize(n) for n in G.nodes()}
    for can_node in unique_canonical_nodes:
        vec = contig_to_freq_vector(SEQs[can_node])
        node_vec_dict[can_node] = vec
    return node_vec_dict

def get_node_scores(scores_file,G):
    """ Write the plasmid scores into each node in the graph
    """
    scores = {}
    unique_canonical_nodes = {canonicalize(n) for n in G.nodes()}
    with open(scores_file) as f:
        for line in f:
            split = line.strip().split()
            if split[0] in unique_canonical_nodes: 
                scores[split[0]] = float(split[1])

    return scores

def get_gene_nodes(genes_file,G):
    """ Annotate each node in the graph whether it has plasmid gene in it
    """
    gene_nodes = set()
    with open(genes_file) as f:
        for line in f:
            gene_nodes.add(line.strip())

    return gene_nodes

def rc_seq(dna):
    rev = reversed(dna)
    return "".join([complements[i] for i in rev])

def get_num_from_spades_name(name):
    name_parts = name.split("_")
    contig_num = name_parts[1]
    return contig_num

def get_length_from_spades_name(name):
    name_parts = name.split("_")
    contig_length = name_parts[3]
    return int(contig_length)

def get_cov_from_spades_name(name):
    name_parts = name.split("_")
    cov = name_parts[5]
    if cov[-1]=="'": cov=cov[:-1]
    return float(cov)

def get_fastg_digraph(fastg_name):
    """ scans through fastg headers as an adjacency list
        builds and returns a nx directed graph using adjacencies
        note: no connections are created between each node and its
        rc node - we need to take care to maintain these
    """
    lines = []
    fp = open(fastg_name, 'r')
    for name,seq,qual in readfq(fp):
        name = re.sub('[:,]'," ", name[:-1])
        lines.append(name)
    G = nx.DiGraph()
    return nx.parse_adjlist(lines, create_using=G)

def get_fastg_seqs_dict(fastg_name, G):
    """ returns a dictionary of sequences in graph
        where node names are keys and sequence strings
        are values; useful for saving memory when G
        is a subgraph (e.g., a component)
    """
    fp = open(fastg_name, 'r')
    seqs = {}
    id_to_fullname = {}
    for name,seq,qual in readfq(fp):
        name_parts = re.sub('[:,]'," ", name[:-1]).split()
        node = name_parts[0]
        seqs[node] = seq
        node_id = node.split("_")[1]
        id_to_fullname[node_id] = node[:-1] if node.endswith("'") else node
    return seqs, id_to_fullname

def rc_node(node):
    """ gets reverse complement
        spades node label
    """
    if node[-1] == "'": return node[:-1]
    else: return node + "'"

def get_cov_from_spades_name_and_graph(name,G):
    if name not in G:
        return 0
    if 'cov' in G.nodes[name]:
        return G.nodes[name]['cov']
    else:
        return get_cov_from_spades_name(name)

def update_node_coverage(G, node, new_cov):
    """ changes coverage value stored in 'cov'
        field on both F and R version of a node
        if new_cov is 0, both node versions are removed
    """
    if node not in G.nodes(): # nothing to be done, perhaps already removed
        return
    if new_cov <1e-9:
        G.remove_node(node)
        logger.info(f"remove {node}")
        if rc_node(node) in G.nodes():
            G.remove_node(rc_node(node))
    else:
        G.add_node(node, cov=new_cov)
        if rc_node(node) in G.nodes():
          G.add_node(rc_node(node), cov=new_cov)

def get_spades_base_mass(G, name):
    length = get_length_from_spades_name(name)
    coverage = get_cov_from_spades_name_and_graph(name,G)
    if coverage <= 1e-9: coverage = 1.0/float(length) # ensure no division by zero, consider more principled way to do this
    return length * coverage

def get_seq_from_path(path, seqs, max_k_val=77, cycle=True):
    """ retrieves sequence from a path;
        instead of specifying cycles by having the first and
        last node be equal, the user must decide if a cycle is
        the intent to avoid redundant k-mers at the ends
    """
    start = seqs[path[0]]
    if len(path)==1:
        if cycle:
            return start[max_k_val:]
        else:
            return start
    else:
        seq = ''
        for p in path:
            seq += seqs[p][max_k_val:]
        if cycle: return seq
        else: return start[:max_k_val] + seq

def get_wgtd_path_coverage_CV(path, G, seqs, max_k_val=77):
    if len(path)< 2: return 0
    mean, std = get_path_mean_std(path, G, seqs, max_k_val,discount=True)
    if mean<=0: return 0
    return std/mean

def get_node_cnts_hist(path):
    d = {}
    for p in path:
        # always count based on positive node
        pos_name = p if (p[-1]!="'") else p[:-1]
        d[pos_name] = d.get(pos_name,0) + 1
    return d

def get_discounted_node_cov(node,path,G):
    """ Return the coverage of the node, discounted by the coverage of neighbouring
        nodes not in the path
    """
    pred_covs = [(get_cov_from_spades_name_and_graph(p,G),p) for p in G.predecessors(node)]
    succ_covs = [(get_cov_from_spades_name_and_graph(s,G),s) for s in G.successors(node)]

    non_path_cov = sum([p[0] for p in pred_covs if p[1] not in path]) + sum([s[0] for s in succ_covs if s[1] not in path])
    in_path_cov = sum([p[0] for p in pred_covs if p[1] in path]) + sum([s[0] for s in succ_covs if s[1] in path])
    node_cov = get_cov_from_spades_name_and_graph(node,G)
    node_cov *= in_path_cov/(non_path_cov + in_path_cov)
    ###################### A possible alternative would be to discount in this way for both the in- and out-neighbours
    ###################### and then average the two discounted
    return node_cov

def get_discounted_node_cov_optimized(node, path, G):
    """ 
    Return the coverage of the node, discounted separately by the coverage 
    of predecessor and successor nodes not in the path, then averaged.
    """
    node_cov = get_cov_from_spades_name_and_graph(node, G)

    # 1. 前驱（Pred）分析
    pred_covs = [(get_cov_from_spades_name_and_graph(p, G), p) for p in G.predecessors(node)]
    
    # 路径中的前驱覆盖度
    in_path_pred_cov = sum([p[0] for p in pred_covs if p[1] in path])
    # 所有前驱覆盖度之和
    all_pred_cov = sum([p[0] for p in pred_covs])
    
    # 计算前驱折扣覆盖度
    if all_pred_cov > 0:
        pred_discount_ratio = in_path_pred_cov / all_pred_cov
        disc_cov_pred = node_cov * pred_discount_ratio
    else:
        # 如果没有前驱，则认为在输入端是完全忠诚的
        disc_cov_pred = node_cov

    # 2. 后继（Succ）分析
    succ_covs = [(get_cov_from_spades_name_and_graph(s, G), s) for s in G.successors(node)]
    
    # 路径中的后继覆盖度
    in_path_succ_cov = sum([s[0] for s in succ_covs if s[1] in path])
    # 所有后继覆盖度之和
    all_succ_cov = sum([s[0] for s in succ_covs])

    # 计算后继折扣覆盖度
    if all_succ_cov > 0:
        succ_discount_ratio = in_path_succ_cov / all_succ_cov
        disc_cov_succ = node_cov * succ_discount_ratio
    else:
        # 如果没有后继，则认为在输出端是完全忠诚的
        disc_cov_succ = node_cov

    # 3. 最终平均
    # 如果节点是路径的起点或终点，它可能只有前驱或后继
    return (disc_cov_pred + disc_cov_succ) / 2.0

def get_path_covs(path, G, discount=False):
    cnts = {}
    if discount:
        # 1. 先统计节点出现次数
        cnts = get_node_cnts_hist(path)
        
        # 2. 使用缓存避免重复计算昂贵的锚点搜索
        # 即使路径很长，唯一节点数通常不多
        node_cov_cache = {} 
        covs = []
        
        for n in path:
            # 规范化节点名用于缓存键值（正负链共享同一个物理覆盖度）
            pos_name = n if (n[-1] != "'") else n[:-1]
            
            if pos_name not in node_cov_cache:
                # 只计算一次
                raw_anchored_cov = calc_anchored_cov(n, path, G)
                
                # 在这里直接做除法，存入缓存
                if cnts[pos_name] > 1:
                    node_cov_cache[pos_name] = raw_anchored_cov / cnts[pos_name]
                else:
                    node_cov_cache[pos_name] = raw_anchored_cov
            
            covs.append(node_cov_cache[pos_name])
            
    else:
        covs = [get_cov_from_spades_name_and_graph(n, G) for n in path]

    return covs

def get_path_mean_std(path, G, seqs, max_k_val=77,discount=True):
    covs = np.array(get_path_covs(path,G,discount))
    wgts = np.array([(get_length_from_spades_name(n)-max_k_val) for n in path])
    tot_len = get_total_len_from_path(path,max_k_val=max_k_val, cycle=True)
    if tot_len<=0: return (0,0)
    wgts = np.multiply(wgts, 1./tot_len)
    mean = np.average(covs, weights = wgts)
    std = np.sqrt(np.dot(wgts,(covs-mean)**2))
    return (mean,std)

def update_path_coverage_vals(path, G, seqs, max_k_val=77,):
    mean, _ = get_path_mean_std(path, G, seqs, max_k_val) ## NOTE: CAN WE STILL GUARANTEE CONVERGENCE WHEN DISCOUNTING COVERAGE ??!
    covs = get_path_covs(path,G)
    discounted_covs = get_path_covs(path,G,True)
    # 如果 repeat的discount_cov过低会拉低整体的mean,此时剥离一个cycle时，不一定能去掉(打断)一个环，不能保证收敛，
    # 后续寻找cycle时，这些cov几近为0但score很高的node会提供错误的信息，误导cycle的寻找
    new_covs = []
    mean = max(mean, min(covs)) if covs else mean
    # 对于多次出现的节点，理应将其折算的部分加进去
    # 如 A(10), B(24), C(10), B'(24)
    # 计算mean时会对相同节点进行折算：
    # ==>covs: A(10), B(12), C(10), B'(12)   mean = 11
    # ==> new_covs: A(0), B(13), C(0), B'(13) 
    # 但是new_covs是直接应用在node上的，即节点B的cov直接被赋予13，但是考虑B这个点，原本为24，删除了两次11，结果理应为2
    # 所以这里需要将new_covs乘以对应节点出现的次数，本质上是因为B和B'在组装图中代表一个顶点
    cnts = get_node_cnts_hist(path)


    for i in range(len(path)):
        p = path[i]
        pos_name = p if (p[-1]!="'") else p[:-1]
        if cnts[pos_name] > 1:
            new_covs.append(covs[i]- mean*cnts[pos_name])
        else:
            new_covs.append(covs[i]-mean)
    
    logger.info("Path: %s Mean: %s Covs: %s DisCountedCovs: %s" % (str(path),str(mean),str(covs),str(discounted_covs) ))
    logger.info("NewCovs: %s" % (str(new_covs)))
    for i in range(len(path)):
        if new_covs[i] > 1e-9:
            update_node_coverage(G,path[i],new_covs[i])
        else:
            update_node_coverage(G,path[i],0)
    return new_covs

def update_path_with_covs(path, G, covs):
    for i in range(len(path)):
        if covs[i] > 1e-9:
            update_node_coverage(G,path[i],covs[i])
        else:
            update_node_coverage(G,path[i],0)

def get_total_path_mass(path,G):
    return sum([get_length_from_spades_name(p) * \
        get_cov_from_spades_name_and_graph(p,G) for p in path])

def get_long_self_loops(G,node_score_dict,node_gene_set, min_length, seqs, bamfile, use_scores=True, use_genes=True, max_k_val=77, score_thresh=0.9, mate_thresh = 0.1):
    """ returns set of self loop nodes paths that are longer
        than min length and satisfy mate pair requirements;
        removes those and short self loops from G
    """
    potential_plasmids = set([])
    to_remove = []

    for nd in list(nx.nodes_with_selfloops(G)):
        if (rc_node(nd),) in potential_plasmids: continue
        nd_path = (nd,)
        path_len = len(get_seq_from_path(nd_path, seqs, max_k_val))

        # check whether it is isolated or connected to other nodes:
        isolated_loop = False
        if G.in_degree(nd) == 1 and G.out_degree(nd)== 1:
            isolated_loop = True
        if isolated_loop:
            if path_len < min_length:
                to_remove.append(nd)
                continue

            # take nodes that have plasmid genes or very high plasmid scores
            if use_scores and use_genes:
                logger.info("SLS: %f" % PARAMS.SELF_LOOP_SCORE_THRESH)
                if node_score_dict[canonicalize(nd)] > PARAMS.SELF_LOOP_SCORE_THRESH or nd in node_gene_set:
                    potential_plasmids.add(nd_path)
                    logger.info("Added path: %s - high scoring long self-loop" % nd)
                    to_remove.append(nd)
                    continue

            off_node_mate_count, on_node_mate_count = count_selfloop_mates(nd,bamfile)
            if float(off_node_mate_count) > PARAMS.SELF_LOOP_MATE_THRESH*float(on_node_mate_count):
                logger.info('Self loop %s has %2f percent off-node mate-pairs. Removing' % (nd,PARAMS.SELF_LOOP_MATE_THRESH))
                to_remove.append(nd)
            else:
                potential_plasmids.add(nd_path)
                logger.info("Added path: %s  - long self loop" % nd)
                to_remove.append(nd)
        else: # non-isolated loop
            if path_len < min_length: continue

            off_node_mate_count, on_node_mate_count = count_selfloop_mates(nd,bamfile)
            if float(off_node_mate_count) > PARAMS.SELF_LOOP_MATE_THRESH*float(on_node_mate_count):  # TODO: could be different than for isolated loop
                                                                                    # Maybe - func of node length (and read length, insert size???)
                logger.info('Self loop %s has %2f percent off-node mate-pairs.' % (nd,PARAMS.SELF_LOOP_MATE_THRESH))
            else:
                potential_plasmids.add(nd_path)
                logger.info("Added path: %s  - long self loop" % nd)
                to_remove.append(nd)

    for nd in to_remove:
        update_node_coverage(G, nd, 0)
    logger.info("Removing %d self-loop nodes" % len(to_remove))
    return potential_plasmids

def remove_hi_confidence_chromosome(G,node_to_contig,score_dict):
    """ Remove the long nodes that are predicted to likely be chromosomal
        Retain nodes in potential plasmid contigs
    """
    to_remove = []
    for nd in G.nodes():
        if get_length_from_spades_name(nd) > PARAMS.CHROMOSOME_LEN_THRESH and \
            score_dict[canonicalize(nd)] < PARAMS.CHROMOSOME_SCORE_THRESH and \
            nd not in node_to_contig:
            to_remove.append(nd)

            to_remove.append(rc_node(nd))
    G.remove_nodes_from(to_remove)
    # logger.info(f"to remove: {str(to_remove)}")
    logger.info("Removed %d long, likely chromosomal nodes" % len(set(to_remove)))

def get_hi_conf_plasmids(G,score_dict):
    """ Return a list of nodes that are likely plasmids
    """

    hi_conf_plasmids = [nd for nd in G.nodes() if (get_length_from_spades_name(nd) > PARAMS.PLASMID_LEN_THRESH and \
                        score_dict[canonicalize(nd)] > PARAMS.PLASMID_SCORE_THRESH)]

    logger.info("Found %d long, likely plasmid nodes" % len(hi_conf_plasmids))
    return hi_conf_plasmids

def get_plasmid_gene_nodes(G,gene_set):
    """ Return list of nodes annotated as having a plasmid gene
    """
    plasmid_gene_nodes = [nd for nd in gene_set if nd in G.nodes()]
    logger.info("Found %d nodes with plasmid genes" % len(plasmid_gene_nodes))
    return plasmid_gene_nodes

def get_unoriented_sorted_str(path):
    """ creates unique, orientation-oblivious string representation of path,
        used to make sure node covered whenever rc of node is;
        lets us avoid issue of rc of node having different weight than node
    """
    all_rc_path = []
    for p in path:
        if p[-1] != "'": p = p+"'"
        all_rc_path.append(p)
    return "".join(sorted(all_rc_path))


def estimate_insert_size_distribution(bamfile, max_samples=10000):
    """
    流式估计 insert size 的均值和标准差，最多采样 max_samples 条记录。
    """
    n = 0
    mean = 0.0
    M2 = 0.0  # 用于计算方差

    for hit in bamfile:
        if n >= max_samples:
            break
        if (hit.is_proper_pair 
            and not hit.is_unmapped 
            and not hit.mate_is_unmapped):
            tlen = abs(hit.template_length)
            if tlen > 0:
                n += 1
                delta = tlen - mean
                mean += delta / n
                delta2 = tlen - mean
                M2 += delta * delta2

    if n < 2:
        return 300, 50  # 默认值

    variance = M2 / (n - 1)  # 无偏估计（sample std）
    std = math.sqrt(variance)
    return mean, std

def get_physical_position(contig, read_pos, read_strand, contig_dir):

    """
    返回 read 在物理 DNA 上的实际位置（相对于 5' → 3'）
    :param contig: contig 名称
    :param read_pos: read 的 mapping 位置（0-based）
    :param read_strand: read 的 mapping 方向 '+' or '-'
    :param contig_dir: contig 在路径中的方向 '+' or '-'
    :return: 物理位置（0-based）
    """
    length = get_length_from_spades_name(contig)

    # read 在 contig 上的“原始”位置（正向 contig 的视角）
    if read_strand == '+':
        raw_pos = read_pos
    else:
        raw_pos = length - read_pos - 1  # reverse read 的原始位置

    # contig 被反向使用时,位置也要翻转
    if contig_dir == '-':
        return length - raw_pos - 1
    else:
        return raw_pos 

def build_pe_support_dict(pe_contigs_path_dict):
    """
    统计每条有向边 (u, v) 被 PE 证据路径支持的次数
    
    :param pe_contigs_path_dict: get_pe_support_evidence() 的返回值
    :return: defaultdict(int), key: (u, v) -> total_score
    """
    pe_support_dict = defaultdict(int)
    total_support = 0

    for direction_dict in pe_contigs_path_dict:  # 遍历正向 [0] 和反向 [1]
        for u, path_list in direction_dict.items():
            for (rest_path, label, score, _) in path_list:
                # 累加整条路径的 score 到总支持数
                total_support += score

                # === 构建路径中的所有相邻边 ===
                path_edges = []

                # 第一条边: u -> rest_path[0]
                if len(rest_path) > 0:
                    path_edges.append((u, rest_path[0]))

                # 后续边: rest_path[i] -> rest_path[i+1]
                for i in range(len(rest_path) - 1):
                    path_edges.append((rest_path[i], rest_path[i + 1]))
                # 累加每条边的支持数
                for edge in path_edges:
                    pe_support_dict[edge] += score  

    logger.info(f"Total PE support count (by path): {total_support}")
    return pe_support_dict

def find_all_paths_with_mid_length_cutoff(
    G, start, target, cutoff, node_length_cache,
    max_hops=10, max_paths=100
):
    if start not in G or target not in G:
        return []
    if start == target:
        return [[start]]

    all_paths = []
    queue = deque()
    queue.append((start, [start], 0, 0))  # (node, path, mid_len, hops)

    while queue and len(all_paths) < max_paths:
        current, path, mid_len, hops = queue.popleft()

        if current == target:
            if mid_len <= cutoff:
                all_paths.append(path)
            continue

        if hops >= max_hops or mid_len > cutoff:
            continue

        for neighbor in G.successors(current):
            if neighbor in path:  # simple path only
                continue

            add_len = 0 if neighbor == target else node_length_cache.get(neighbor, 0)
            new_mid_len = mid_len + add_len
            new_hops = hops + 1

            if new_mid_len > cutoff and neighbor != target:
                continue

            queue.append((neighbor, path + [neighbor], new_mid_len, new_hops))

    return all_paths


# --- Main function ---
def get_pe_support_evidence(G, bamfile, insert_mean, insert_std, max_k=77, max_hops=10, max_paths_per_pair=50):
    """
    Extract PE support evidence from BAM file.
    Returns:
        pe_contigs_path_dict[0]: forward paths
        pe_contigs_path_dict[1]: reverse paths
        valid_mate_pairs: dict of mutual supports
    """
    logger.info(f"k = {max_k}, insert mean = {insert_mean}, std = {insert_std}")

    # === Step 0: Precompute global info ===
    # Node length cache
    node_length_cache = {}
    for node in G.nodes():
        base = canonicalize(node)
        node_length_cache[node] = get_length_from_spades_name(base)

    # Component map for fast connectivity check
    comp_id_map = {}
    for cid, comp in enumerate(nx.connected_components(G.to_undirected())):
        for node in comp:
            comp_id_map[node] = cid

    # Result containers
    pe_contigs_path_dict = [defaultdict(list), defaultdict(list)]
    valid_mate_pairs = defaultdict(set)

    # === Step 1: Stream and group PE pairs by (u, v) ===
    logger.info("Streaming and grouping valid PE reads...")
    pe_grouped = defaultdict(list)
    seen_read_ids = set()

    for read in bamfile.fetch():
        if not (read.is_paired and not read.is_unmapped and not read.mate_is_unmapped):
            continue
        if read.is_read2:
            continue

        qname = read.query_name
        if qname in seen_read_ids:
            continue
        seen_read_ids.add(qname)

        u_contig = read.reference_name
        v_contig = bamfile.getrname(read.next_reference_id)
        if u_contig is None or v_contig is None or u_contig == v_contig:
            continue

        # Strand & position
        u_strand = '-' if read.is_reverse else '+'
        v_strand = '-' if read.mate_is_reverse else '+'

        u_contig_strand = '-' if u_contig.endswith("'") else '+'
        v_contig_strand = '-' if v_contig.endswith("'") else '+'

        u_pos = get_physical_position(u_contig, read.reference_start, u_strand, u_contig_strand)
        v_pos = get_physical_position(v_contig, read.next_reference_start, v_strand, v_contig_strand)

        pe_grouped[(u_contig, v_contig)].append({
            'u': u_contig,
            'v': v_contig,
            'u_pos': u_pos,
            'v_pos': v_pos,
        })

    logger.info(f"Grouped into {len(pe_grouped)} (u,v) pairs.")

    # === Step 2: Process each (u, v) group ===
    total_paths_found = 0
    cutoff = insert_mean + 2 * insert_std

    for (u_contig, v_contig), pairs in pe_grouped.items():
        u_comp = comp_id_map.get(u_contig)
        v_comp = comp_id_map.get(v_contig)
        if u_comp is None or v_comp is None or u_comp != v_comp:
            continue  # skip if not in same component

        for is_rc in [False, True]:
            s = rc_node(u_contig) if is_rc else u_contig
            t = rc_node(v_contig) if is_rc else v_contig

            if s not in G or t not in G:
                continue

            # Fast path: direct edge?
            if G.has_edge(s, t):
                paths = [[s, t]]
            else:
                paths = find_all_paths_with_mid_length_cutoff(
                    G, s, t, cutoff, node_length_cache,
                    max_hops=max_hops, max_paths=max_paths_per_pair
                )

            if not paths:
                continue

            v_len = node_length_cache.get(v_contig, 0)

            for path in paths:
                total_len = sum(
                    node_length_cache.get(node, 0) - (max_k - 1) if i > 0 else node_length_cache.get(node, 0)
                    for i, node in enumerate(path)
                )
                # Simpler: total physical length of path
                # If you have a helper, use: total_len = get_total_len_from_path(path, max_k, cycle=False)

                # Count valid pairs based on insert size estimate
                valid_count = 0
                for p in pairs:
                    # insert_size ≈ u_pos + (path length) + (v_len - v_pos)
                    insert_est = p['u_pos'] + total_len + (v_len - p['v_pos'])
                    if abs(insert_est - insert_mean) <= 2 * insert_std:
                        valid_count += 1

                if valid_count < 2:
                    continue

                # Store forward path
                pe_contigs_path_dict[0][path[0]].append((
                    tuple(path[1:]), "mate", valid_count, None
                ))

                # Store reverse path
                rev_path = list(reversed(path))
                pe_contigs_path_dict[1][rev_path[0]].append((
                    tuple(rev_path[1:]), "Rmate", valid_count, None
                ))

                valid_mate_pairs[s].add(t)
                valid_mate_pairs[t].add(s)
                total_paths_found += 1

    logger.info(f"Found {total_paths_found} PE support path instances.")
    return pe_contigs_path_dict, dict(valid_mate_pairs)

def get_weighted_cov(G, in_nodes:list, out_nodes:list, t: str, max_k_val: int):

    def get_contig_path_mean(path, G, max_k_val=77,discount=True):

        if not path:
            return 0.0
            
        covs = []
        lengths = []
        
        for node in path:
            # 获取原始覆盖度（不需要锚点折算，因为我们看的是节点本身的丰度）
            c = get_cov_from_spades_name_and_graph(node, G)
            # 获取有效长度 (len - k + 1)
            l = get_length_from_spades_name(node) - max_k_val + 1
            if l < 1: l = 1
            
            covs.append(c)
            lengths.append(l)
            
        
        covs = np.array(covs)
        lengths = np.array(lengths)
        
      
        sorted_indices = np.argsort(covs)
        sorted_covs = covs[sorted_indices]
        sorted_weights = lengths[sorted_indices]
        

        cumsum_weights = np.cumsum(sorted_weights)
        total_weight = cumsum_weights[-1]
        

        median_idx = np.searchsorted(cumsum_weights, total_weight / 2.0)
        
        return sorted_covs[median_idx]
    
    if t not in ["in", "out"]:
        raise ValueError("type must be 'in' or 'out'")
    
    def to_node_list(x):
        return [x] if isinstance(x, str) else list(x)
    
    in_nodes = to_node_list(in_nodes)
    out_nodes = to_node_list(out_nodes)

    cov_in = get_cov_from_spades_name_and_graph(in_nodes[0], G) if len(in_nodes) == 1 else get_contig_path_mean(in_nodes,G,max_k_val)
    cov_out= get_cov_from_spades_name_and_graph(out_nodes[0], G) if len(out_nodes) == 1 else get_contig_path_mean(out_nodes,G,max_k_val)


    if t == "in":
        neighbors = G.successors(in_nodes[-1])
    else:  # out
        neighbors = G.predecessors(out_nodes[0])
    
    # 收集边权重
    weights = []
    
    for neighbor in neighbors:
        weight = get_cov_from_spades_name_and_graph(neighbor,G)
        weights.append(weight)
    
    
    total_weight = sum(weights)
    
    # 如果邻接顶点属于contig path,则补偿总权值（因为cov_out已经包含了out_nodes[0]）
    if t == "in" and len(out_nodes) > 1:
        total_weight += cov_out - get_cov_from_spades_name_and_graph(out_nodes[0],G)
    if t == "out" and len(in_nodes) > 1:
        total_weight += cov_in - get_cov_from_spades_name_and_graph(in_nodes[-1],G)

    if total_weight == 0:
        raise ValueError("Total weight is zero, cannot compute coverage")
    


    return cov_out * cov_in / total_weight

def get_supports(G: nx.DiGraph, in_nodes: list, out_nodes: list):
    """
    计算 in_nodes 中所有节点到 out_nodes 中所有节点的 PE 总支持数。
    从图 G 的属性中获取 pe_support_dict。
    """
    # 1. 从图属性中获取 PE 支持字典
    # 使用 .get() 避免在图 G 中没有该属性时出错
    pe_support_dict = G.graph.get('pe_support_data', {}) 
    
    # 如果字典不存在，直接返回 0
    if not pe_support_dict:
        return 0

    # support = 0
    # # for u in in_nodes:
    # #     for v in out_nodes:
    # #         support += pe_support_dict.get((u, v), 0)
    
    return pe_support_dict.get((in_nodes[-1], out_nodes[0]), 0)


def get_edge_cost(G, in_nodes:list, out_nodes:list, in_score:int, out_score:int, in_vec, out_vec,support, max_k: int):
   
    
    in_cov = get_weighted_cov(G, in_nodes, out_nodes, t='in', max_k_val = max_k)
    out_cov = get_weighted_cov(G, in_nodes, out_nodes, t='out', max_k_val = max_k)
    

    cost = (1 - (np.sqrt(in_score * out_score))) + \
           (1 - cosine_similarity(in_vec.reshape(1, -1), out_vec.reshape(1, -1))[0][0]) + \
           abs(in_cov - out_cov) / max(in_cov, out_cov) + \
           (1 / (1 + support)) 
           
    return cost

def get_shortest(args_array):
    """ Worker function for getting shortest path to each node in parallel
    """
    
    node, path_dict,node_score_dict,node_vec_dict,node_support_dict,SEQS,G = args_array
    shortest_score = float("inf")
    path = None
    use_contig = False
    paths_list = []

    for pred in G.predecessors(node):
        try:
            path_len,shortest_path,use = dijkstra_path(G,path_dict,node_score_dict,node_vec_dict,node_support_dict,SEQS,source=node,target=pred, weight='cost',bidirectional=True)
            
            # logger.info(f"shortest_path:  {str(folded_path)}")
            if path_len < shortest_score:
                path = shortest_path
                shortest_score = path_len
                use_contig = use
        except nx.exception.NetworkXNoPath:
            continue
    # logger.info(f"add shortest path: {path}")
    if path is not None: paths_list.append((path,shortest_score,use_contig))
    return paths_list
    # done


def enum_high_mass_shortest_paths(G,path_dict,node_gene_set,node_score_dict,node_vec_dict,node_support_dict,SEQS, use_scores=False, use_genes=False, seen_paths=None, max_k=77,batch_size=200,num_procs=8):
    """ given component subgraph, returns list of paths that
        - is non-redundant (includes) no repeats of same cycle
        - includes shortest paths starting at each node n (assigning
        node weights to be 1/(length * coverage)) to each of
        its predecessors, and returning to n
    """
    current_nodes = set(G.nodes())  # 当前连通分量的所有节点

    def get_canonical_path(path):
        rev_path = tuple(rc_node(node) for node in reversed(path))
        forward_exists = all(node in current_nodes for node in path)
        reverse_exists = all(node in current_nodes for node in rev_path)
        
        if forward_exists and reverse_exists:
            return path if path <= rev_path else rev_path
        elif forward_exists:
            return path
        elif reverse_exists:
            return rev_path
        

    if seen_paths == None:
        seen_paths = []
    unq_sorted_paths = set([])
    # in case orientation obliv. sorted path strings passed in
    for p in seen_paths:
        unq_sorted_paths.add(p)
    paths = []

    # logger.info("Getting edge weights")

    # use add_edge to assign edge weights to be 1/mass of starting node
    # TODO: only calculate these if they haven't been/need to be updated
    for e in G.edges():
        # if use_genes and G.nodes[e[1]]['gene'] == True:
        #     G.add_edge(e[0], e[1], cost = 0.0)
        # elif use_scores==True:
           G.add_edge(e[0], e[1], cost = get_edge_cost(G, e[0], e[1], node_score_dict[canonicalize(e[0])], node_score_dict[canonicalize(e[1])],\
                                                       node_vec_dict[canonicalize(e[0])], node_vec_dict[canonicalize(e[1])],node_support_dict.get((e[0],e[1]),0),max_k=max_k))
        # else:
            # G.add_edge(e[0], e[1], cost = (1./get_spades_base_mass(G, e[1])))
    valid_path_starts = set(path_dict[0].keys())
    logger.info("Getting shortest paths")
    nodes = [n for n in current_nodes if (get_length_from_spades_name(n) >= 1000 \
             or n in valid_path_starts\
             or n in node_gene_set)
            ]
    remain_nodes = len(nodes)
    print(str(remain_nodes) + " nodes remain in component")
    logger.info("Remaining nodes: %d" % (remain_nodes))
    batches = [nodes[i:i + batch_size] for i in range(0, remain_nodes, batch_size)]
    

    use_parallel = True if num_procs >1 else False 

    if use_parallel and len(batches)>1:
        local_pool = mp.Pool(
            processes=num_procs, 
            initializer=init_worker, 
            initargs=(path_dict, node_score_dict, node_vec_dict, node_support_dict, SEQS, G)
        )
        logger.info(f"use multi precessor with {len(batches)} batches)")
        try:
            results = local_pool.map(get_shortest_batch, batches)
        finally:
            local_pool.close()
            local_pool.join()
    else:
        logger.info(f"use single precessor")
        init_worker(path_dict, node_score_dict, node_vec_dict, node_support_dict, SEQS, G)
        results = [get_shortest_batch(batch) for batch in batches]

    # --- 合并结果 ---
    paths_list = []
    for res in results:
        if res:
            paths_list.extend(res)

    paths_list.sort(key=sort_key)
    # --- 去重（结合 seen_paths）---

    path_to_seed_dict = defaultdict(list)

    final_paths = []
    for path_tuple, weight, use_contig, seed in paths_list:
        # 生成无向路径字符串用于去重
        path_tuple = get_canonical_path(path_tuple)
        unoriented_str = get_unoriented_sorted_str(path_tuple)
        path_to_seed_dict[path_tuple].append(seed)
        if unoriented_str not in unq_sorted_paths:
            unq_sorted_paths.add(unoriented_str)
            final_paths.append((path_tuple, weight, use_contig))


    return final_paths, path_to_seed_dict


def get_high_mass_shortest_path(node,G,path_dict,proxy_contig_dict,node_score_dict,node_vec_dict,node_support_dict,SEQS,use_scores,use_genes,max_k):
    """ Return the shortest circular path back to node
    """
    # rc_node_hit used to determine whether use node or rc_node(node) to find plasmid gene hit cycle

    # TODO: potentially add check for unique paths so that don't check same cycle
    # twice if there are two potential plasmid nodes in it

    for e in G.edges():
        # if use_genes and G.nodes[e[1]]['gene'] == True:
        #     G.add_edge(e[0], e[1], cost = 0.0)
        # elif use_scores == True:
            G.add_edge(e[0], e[1], cost = get_edge_cost(G, e[0], e[1], node_score_dict[canonicalize(e[0])], node_score_dict[canonicalize(e[1])],\
                                                        node_vec_dict[canonicalize(e[0])], node_vec_dict[canonicalize(e[1])],node_support_dict.get((e[0],e[1]),0),max_k=max_k))
        # else:
            # G.add_edge(e[0], e[1], cost = (1./get_spades_base_mass(G, e[1])))

    shortest_score = float("inf")
    path = None
    shortest_path = None
    use_contig = None
    if proxy_contig_dict:
        if node in proxy_contig_dict.keys():
            for path_info in proxy_contig_dict[node]:
                path, name, score, pre,freq_vec = path_info
                if pre is not None:
                    proxy = pre[0]
                    if proxy not in G.nodes(): continue
                    for pred in G.predecessors(proxy):
                        try:
                            length, path, through_contig = dijkstra_path(G,path_dict,node_score_dict,node_vec_dict,node_support_dict,SEQS,source=proxy,target=pred, weight='cost')
                            if length < shortest_score:
                                shortest_path = tuple(path)
                                shortest_score = length
                                use_contig = through_contig
                            # logger.info(f"special path:  {str(folded_path)}")
                        except nx.exception.NetworkXNoPath:
                            continue
    for pred in G.predecessors(node):
        try:
            length, path, through_contig = dijkstra_path(G,path_dict,node_score_dict,node_vec_dict,node_support_dict,SEQS,source=node,target=pred, weight='cost')
            if length < shortest_score:
                shortest_path = tuple(path)
                shortest_score = length
                use_contig = through_contig
            # logger.info(f"special path:  {str(folded_path)}")
        except nx.exception.NetworkXNoPath:
            continue

    return shortest_score,shortest_path,use_contig

def get_non_repeat_nodes(G, path):
    """ returns a list of all non-repeat (in degree and out-degree
        == 1) nodes in a path; if there are no such nodes,
        returns an empty list
        NB: G input should be whole graph, not specific SCC, to avoid
        disregarding isolated nodes
    """
    sing_nodes = []
    for nd in path:
        if G.out_degree(nd)==1 and G.in_degree(nd)==1:
            sing_nodes.append(nd)
    return sing_nodes

def get_spades_type_name(count, path, seqs, max_k_val, G, cov=None):
    path_len = len(get_seq_from_path(path,seqs,max_k_val))
    if cov==None:
        cov = get_total_path_mass(path,G)/float(path_len)
    info = ["RNODE", str(count+1), "length", str(path_len),
     "cov", '%.5f' % (cov)]
    return "_".join(info)

def count_selfloop_mates(node,bamfile):
    """ Counts the number of off-node and on-node mate pairs of
        a self-loop node
    """
    off_node_count = 0
    on_node_count = 0
    if node[-1] == "'": node = node[:-1]
    try:
        for hit in bamfile.fetch(node):
            nref = bamfile.getrname(hit.next_reference_id)
            if nref != node:
                off_node_count += 1
            else: on_node_count += 1

    except ValueError:
        pass

    return off_node_count, on_node_count

def get_non_path_domain_node(path,valid_mate_pairs):
    """ check all non-repeat nodes only have mates
        mapping to contigs in the cycle
    """
    non_path_dominated_nodes = 0

    for nd in path:
        mate_tigs = valid_mate_pairs.get(nd,set())
      
        num_mates_in_path = sum([1 for x in mate_tigs if (x in path)])
        num_mates_not_in_path = len(mate_tigs)-num_mates_in_path
        if num_mates_in_path < num_mates_not_in_path:
            non_path_dominated_nodes += 1
    return non_path_dominated_nodes

def get_contigs_of_mates(node, bamfile, G):
    """ retrieves set of nodes mapped to by read pairs
        having one mate on node; discards isolated nodes
        because they tend to reflect irrelevant alignments
    """
    
    mate_tigs = set([])
    if node[-1] == "'": node=node[:-1]
    try:
        for hit in bamfile.fetch(node):
            nref = bamfile.getrname(hit.next_reference_id)
            if nref != node:
                mate_tigs.add(nref)

    except ValueError:
        pass

    to_remove = set([])
    for nd in mate_tigs:
        if (G.in_degree(nd)==0 and G.out_degree(nd)==0) or \
        (not G.has_node(nd)):
            to_remove.add(nd)
        # see if nd reachable by node or vice-versa
        # try both flipping to rc and switching source and target
        elif G.has_node(rc_node(node)) and not any([nx.has_path(G, node, nd), nx.has_path(G, rc_node(node),nd),
                nx.has_path(G, nd, node), nx.has_path(G, nd, rc_node(node))]):
            to_remove.add(nd)
        elif not any([nx.has_path(G,node,nd),nx.has_path(G,nd,node)]):
            to_remove.add(nd)
    mate_tigs -= to_remove

    return mate_tigs

def is_good_cyc(path, valid_mate_pairs):
    """ check all non-repeat nodes only have mates
        mapping to contigs in the cycle
    """
    mate_domain_node = 0
    for nd in path:
        mate_tigs = valid_mate_pairs.get(nd, set())
      
        num_mates_in_path = sum([1 for x in mate_tigs if (x in path)])
        num_mates_not_in_path = len(mate_tigs)-num_mates_in_path
        if num_mates_in_path >= num_mates_not_in_path:
            mate_domain_node+=1
    if mate_domain_node < len(path)/2:
        logger.info(f"mate_domain_node: {mate_domain_node}, path_len: {len(path)}")
        return False
    return True
    
def process_component(COMP, G, max_k, min_length, max_CV, SEQS, path_dict,node_to_contig,contigs_path_name_dict,
                      proxy_contig_dict, valid_pairs,node_score_dict,node_gene_set,node_vec_dict,node_support_dict,
                      use_scores=False, use_genes=False, num_procs=1):
    """ run recycler for a single component of the graph
        use multiprocessing to process components in parallel
    """

    # initialize shortest path set considered
    path_count = 0
    seen_unoriented_paths = set([])
    paths_set = set([]) #the set of paths found


    # looking for paths starting from the nodes annotated with plasmid genes
    if use_genes:
        plasmid_gene_nodes = get_plasmid_gene_nodes(COMP,node_gene_set)
        potential_plasmid_mass_tuples = [(get_spades_base_mass(COMP,nd),nd) for nd in plasmid_gene_nodes]
        potential_plasmid_mass_tuples.sort(
            key=lambda n: (n[0], 0 if n[1][-1] == "'" else 1)
        ) # prefer non-rc nodes if masses equal
        while potential_plasmid_mass_tuples: # could be removing other nodes from the list
            top_node = potential_plasmid_mass_tuples.pop() # highest mass node
            top_node_name = top_node[1]
            logger.info(f"plasmid gene hit node: {top_node_name}")
            # 初始化rc_node的相关信息
            rc_length, rc_path,rc_use_contig,rc_path_CV = float("inf"), None, False,float("inf")
            length, path,use_contig = get_high_mass_shortest_path(top_node_name,COMP,path_dict,proxy_contig_dict,node_score_dict,node_vec_dict,node_support_dict,SEQS,use_scores,use_genes,max_k) #######
            path_CV = float("inf")
            if path is not None:
                path_CV = get_wgtd_path_coverage_CV(path,G,SEQS,max_k_val=max_k)
                logger.info(f"plasmid gene path: {path}")
                logger.info(f"CV: {path_CV}, Good: {is_good_cyc(path,valid_pairs)}, weight: {length}, use_contig: {use_contig}")
            if rc_node(top_node_name) in plasmid_gene_nodes:
                rc_length, rc_path,rc_use_contig = get_high_mass_shortest_path(rc_node(top_node_name),COMP,path_dict,proxy_contig_dict,node_score_dict,node_vec_dict,node_support_dict,SEQS,use_scores,use_genes,max_k)
                if rc_path is not None:
                    rc_path_CV = get_wgtd_path_coverage_CV(rc_path,G,SEQS,max_k_val=max_k)
                    logger.info(f"RC plasmid gene path: {rc_path}") 
                    logger.info(f"CV: {rc_path_CV}, Good: {is_good_cyc(rc_path,valid_pairs)}, weight: {rc_length}, use_contig:{rc_use_contig}")
            
            if path is None and rc_path is None: continue
            # 由于正向序列或反向互补序列第一次被添加后就会删除另一个，因此需要通过逻辑判断确定最终的结果：
            # 如果反向互补序列的最短路径更短，并且符合添加条件，则使用反向互补序列，否则使用正向序列
            if path is None:
                path = rc_path
            elif rc_path is not None:
                if meet_criterion(path, G, SEQS,max_k, max_CV,valid_pairs):
                    if meet_criterion(rc_path, G, SEQS,max_k, max_CV,valid_pairs):
                        if rc_use_contig > use_contig or (rc_use_contig == use_contig and rc_path_CV < path_CV):
                            path = rc_path
                elif meet_criterion(rc_path, G, SEQS,max_k, max_CV,valid_pairs):
                    path = rc_path

            if meet_criterion(path, G, SEQS,max_k, max_CV,valid_pairs):
                logger.info("Added plasmid gene path %s" % (str(path)))
                # prevent checking nodes that have been removed
                i = 0
                while i < len(potential_plasmid_mass_tuples):
                    if potential_plasmid_mass_tuples[i][1] in path or \
                        rc_node(potential_plasmid_mass_tuples[i][1]) in path:
                        potential_plasmid_mass_tuples.pop(i)
                    else: i += 1

                seen_unoriented_paths.add(get_unoriented_sorted_str(path))
                # 通过原图G,计算 mean coverage
                before_cov, _ = get_path_mean_std(path, G, SEQS, max_k)
                # 返回path中的节点更新后的coverage
                covs = update_path_coverage_vals(path, G, SEQS, max_k)
                # COMP 利用上述coverage更新自己的coverage
                update_path_with_covs(path, COMP, covs)
                path_count += 1
                paths_set.add((path,before_cov,"high_conf",(path[0],)))

                # 更新contig path information
                path_dict = get_native_path_dict(COMP, path_dict)
                proxy_contig_dict = get_native_proxy_path_dict(COMP, proxy_contig_dict)
            else:
                logger.info("Did not add plasmid gene path: %s" % (str(path)))

                
                
            
        # then look for circular paths that start from hi confidence plasmid nodes
    if use_scores:
        potential_plasmid_nodes = get_hi_conf_plasmids(COMP,node_score_dict)
        potential_plasmid_mass_tuples = [(get_spades_base_mass(COMP,nd),nd) for nd in potential_plasmid_nodes]
        potential_plasmid_mass_tuples.sort(
            key=lambda n: (n[0], 0 if n[1][-1] == "'" else 1)
        ) # prefer non-rc nodes if masses equal
        while potential_plasmid_mass_tuples: # could be removing other nodes from the list
            top_node = potential_plasmid_mass_tuples.pop() # highest mass node
            top_node_name = top_node[1]
            logger.info(f"Hi conf node: {top_node_name}")
            rc_length, rc_path,rc_use_contig,rc_path_CV = float("inf"), None, False,float("inf")
            length, path,use_contig = get_high_mass_shortest_path(top_node_name,COMP,path_dict,proxy_contig_dict,node_score_dict,node_vec_dict,node_support_dict,SEQS,use_scores,use_genes,max_k) #######
            path_CV = float("inf")
            if path is not None:
                path_CV = get_wgtd_path_coverage_CV(path,G,SEQS,max_k_val=max_k)
                logger.info(f"Hi conf path: {path}") 
                logger.info(f"CV: {path_CV}, Good: {is_good_cyc(path,valid_pairs)}, weight: {length}, use_contig: {use_contig}")
                
            if rc_node(top_node_name) in potential_plasmid_nodes:
                rc_length, rc_path,rc_use_contig = get_high_mass_shortest_path(rc_node(top_node_name),COMP,path_dict,proxy_contig_dict,node_score_dict,node_vec_dict,node_support_dict,SEQS,use_scores,use_genes,max_k)
                if rc_path is not None:
                    rc_path_CV = get_wgtd_path_coverage_CV(rc_path,G,SEQS,max_k_val=max_k)
                    logger.info(f"RC Hi conf path: {rc_path}")
                    logger.info(f"CV: {rc_path_CV}, Good: {is_good_cyc(rc_path,valid_pairs)}, weight: {rc_length}, use_contig: {rc_use_contig}")
            
            if path is None and rc_path is None: continue
            # 由于正向序列或反向互补序列第一次被添加后就会删除另一个，因此需要通过逻辑判断确定最终的结果：
            # 如果反向互补序列的最短路径更短，并且符合添加条件，则使用反向互补序列，否则使用正向序列
            if path is None:
                path = rc_path
            elif rc_path is not None and path is not None:
                if meet_criterion(path, G, SEQS,max_k, max_CV,valid_pairs):
                    if meet_criterion(rc_path, G, SEQS,max_k, max_CV,valid_pairs):
                        if rc_use_contig > use_contig or (rc_use_contig == use_contig and rc_path_CV < path_CV):
                            path = rc_path
                elif meet_criterion(rc_path, G, SEQS,max_k, max_CV,valid_pairs):
                    path = rc_path         
                    
            if meet_criterion(path, G, SEQS,max_k, max_CV,valid_pairs):
                logger.info("Added Hi conf path %s" % (str(path)))
                # prevent checking nodes that have been removed
                i = 0
                while i < len(potential_plasmid_mass_tuples):
                    if potential_plasmid_mass_tuples[i][1] in path or \
                        rc_node(potential_plasmid_mass_tuples[i][1]) in path:
                        potential_plasmid_mass_tuples.pop(i)
                    else: i += 1

                seen_unoriented_paths.add(get_unoriented_sorted_str(path))
                before_cov, _ = get_path_mean_std(path, G, SEQS, max_k)
                covs = update_path_coverage_vals(path, G, SEQS, max_k)
                update_path_with_covs(path, COMP, covs)
                path_count += 1
                paths_set.add((path,before_cov,"high_conf",(path[0],)))
                # 更新contig path information
                path_dict = get_native_path_dict(COMP, path_dict)
                proxy_contig_dict = get_native_proxy_path_dict(COMP, proxy_contig_dict)
            else:
                logger.info("Did not add Hi conf path: %s" % (str(path)))


        # 3rd step. Run Recycler algorithm that looks for circular high mass shortest
        # paths and accept them as plasmid predictions if the coverages and mate pairs
        # match the required thresholds
#######################################################################################
#######################################################################################

    paths, seed_to_path = enum_high_mass_shortest_paths(COMP,path_dict,node_gene_set,node_score_dict,node_vec_dict,node_support_dict,SEQS,use_scores,use_genes,seen_unoriented_paths,max_k,num_procs)
    last_path_count = 0
    last_node_count = 0

        # continue as long as you either removed a low mass path
        # from the component or added a new path to final paths
    while(path_count!=last_path_count or\
        len(COMP.nodes())!=last_node_count):

        last_node_count = len(COMP.nodes())
        last_path_count = path_count
        # make tuples of (CV, path)
        path_tuples = []
        path_tuples_with_contig = []
        path_tuples_without_contig = []
        for p,length,use_contig in paths:
            # if len(get_seq_from_path(p, SEQS, max_k_val=max_k)) < min_length:
            if get_total_len_from_path(p,max_k_val=max_k,cycle=True) < min_length:
                seen_unoriented_paths.add(get_unoriented_sorted_str(p))
                # logger.info("Num seen paths: %d" % (len(seen_unoriented_paths)))
                continue
            CV = get_wgtd_path_coverage_CV(p,G,SEQS,max_k_val=max_k)
            if CV > max_CV:continue
            if use_contig:
                path_tuples_with_contig.append((CV, p,length,use_contig))
            else:
                path_tuples_without_contig.append((CV, p,length,use_contig))
            
        if(len(path_tuples_with_contig) > 0):
            path_tuples = path_tuples_with_contig
            logger.info("Using %d paths with contig paths" % (len(path_tuples_with_contig)))
        else:
            path_tuples = path_tuples_without_contig
            logger.info("Using %d paths without contig paths" % (len(path_tuples_without_contig)))
        if(len(path_tuples)==0): break

        # 使用了contig path的优先，低CV作为第二排序依据
        path_tuples.sort(key=lambda item: item[0])

        # logger.info("candidate path:")
        # for cv,p,weight,use_contig in path_tuples:
        #     logger.info(f"\tpath: {p}, CV: {cv}, weight: {weight},use_contig: {use_contig} \n")

        for pt in path_tuples:
            curr_path = pt[1]
            curr_path_CV = pt[0]
            logger.info("Path: %s" % (",".join(curr_path)))
            if get_unoriented_sorted_str(curr_path) not in seen_unoriented_paths:

                ## only report if low CV and matches mate pair info
                if (curr_path_CV <= (max_CV) and \
                    is_good_cyc(curr_path,valid_pairs)):

                    logger.info("Added path %s" % ", ".join(curr_path))
                    logger.info("\tCV: %4f" % curr_path_CV)
                    seen_unoriented_paths.add(get_unoriented_sorted_str(curr_path))
                    before_cov, _ = get_path_mean_std(curr_path, G, SEQS, max_k)
                    covs = update_path_coverage_vals(curr_path, G, SEQS, max_k)
                    update_path_with_covs(curr_path, COMP, covs)
                    path_count += 1
                    paths_set.add((curr_path, before_cov,"remaining",tuple(seed_to_path[curr_path])))

                    # 更新contig path information，这里只需要更新path_dict
                    path_dict = get_native_path_dict(COMP, path_dict)

                    break

                else:
                    logger.info("Did not add path: %s" % (", ".join(curr_path)))
                    logger.info("\tCV: %4f" % curr_path_CV)
                    if curr_path_CV > max_CV:
                        #TODO if there is contig path tuples, we can continue to precess normal tuples(more time-consuming ,more FP and more TP)
                        break # sorted by CV
                    else: # not good mate pairs
                        seen_unoriented_paths.add(get_unoriented_sorted_str(curr_path))

        # recalculate paths on the component
        # print(str(len(COMP.nodes())) + " nodes remain in component")
        # logger.info("Remaining nodes: %d" % (len(COMP.nodes())))
        paths, seed_to_path = enum_high_mass_shortest_paths(COMP,path_dict,node_gene_set,node_score_dict,node_vec_dict,node_support_dict,SEQS, use_scores,use_genes,seen_unoriented_paths,max_k,num_procs)
    # #end while
    st = time.time()
    # path_set_before_merge = paths_set.copy()
    logger.info("start merge...")
    merged_paths_set = merge_cycle(paths_set,SEQS,max_k,node_to_contig,contigs_path_name_dict,valid_pairs)
    logger.info("merge finished:")
    # merged_paths_set = paths_set
    ed = time.time()
    merge_time = ed - st
    return merged_paths_set, merge_time,
    
def merge_cycle(paths_set:set,SEQS,max_k,node_to_contig,contigs_path_name_dict,valid_mate_pairs):
    #记录原有set规模，如果有合并，则规模改变
    prev_len = -1
    #由于此component为强连通分支，如果正向和反向互补序列均在图中，此时选择的cycle已经是两者中更好的，如果反向互补序列不在此连通分支中，则它们理应不连通，也无需考虑
    score_dict = get_score_from_set(paths_set,SEQS)
    while(prev_len != len(paths_set)):
        prev_len = len(paths_set)
        sorted_paths_list = sorted(paths_set,key=lambda x:x[1],reverse=True)
        merged_in_this_iteration = False

        for i in range(len(sorted_paths_list)):
            cur_path,cur_cov, cur_path_type,cur_path_seeds  = sorted_paths_list[i]
            cur_set = set(cur_path)
            cur_score = score_dict[(cur_path,cur_cov)]

            for j in range(i+1, len(sorted_paths_list)):
                other_path, other_cov,other_path_type, other_path_seeds = sorted_paths_list[j]
                original_other_path = other_path
                other_set = set(other_path)
                other_score = score_dict[(other_path,other_cov)]
                intersection = cur_set.intersection(other_set)

                cur_vec = None
                other_vec = None
                similarity = -1

                # 只有当两个路径有交集时才尝试合并
                if len(intersection)>0:
                    # 相似度判断
                    to_merge =  is_similar(cur_cov,other_cov)
                    cur_counter = Counter(cur_path)
                    other_counter = Counter(other_path)
                    
                    diff_cur = cur_counter - other_counter
                    diff_other = other_counter - cur_counter

                    cur_vec = contig_set_to_freq_vector([SEQS[nd] for nd in diff_cur.elements()])
                    other_vec = contig_set_to_freq_vector(SEQS[nd] for nd in diff_other.elements())
                    similarity = cosine_similarity(cur_vec.reshape(1, -1),other_vec.reshape(1, -1))[0][0]
                    merged_path = None
                    logger.info(f"Checking merge: \nCur: {cur_path}  \nOther: {other_path}  \nIntersection: {intersection} \nSimilarity: {similarity}")
                    if to_merge is False:
                        logger.info(f"Not merging based on coverage similarity: cur_cov:{cur_cov} other_cov: {other_cov}")
                        break
                    if similarity < 0.95:
                        logger.info(f"Not merging based on vector similarity:  similarity: {similarity}")
                        continue
                    
                    for nd in sorted(intersection, key=lambda nd: get_cov_from_spades_name(nd), reverse=True):
                        if nd in node_to_contig.keys():
                            merged_path = merge_path_through_contig_path(cur_path,other_path,nd, node_to_contig,contigs_path_name_dict)
                            if merged_path is not None:
                                break

                    if merged_path is not None:
                        merged_score = get_score_from_path(merged_path,SEQS,max_k)
                        
                        logger.info(f"cur: cov: {cur_cov}, score: {score_dict[(cur_path,cur_cov)]}")
                        logger.info(f"other: cov: {other_cov}, score: {score_dict[(other_path,other_cov)]}")
                        merged_cov  = (cur_cov + other_cov)/2
                        logger.info(f"merged_path_through_contig: {merged_path}, cov: {merged_cov}, score: {merged_score}")
                            

                        if merged_score >= 0.5:
                            paths_set.add((merged_path,merged_cov,str(cur_path_type+other_path_type), cur_path_seeds+other_path_seeds))
                            score_dict[(merged_path,merged_cov)] = merged_score
                            #删除原有path
                            paths_set.remove((cur_path,cur_cov,cur_path_type,cur_path_seeds))
                            paths_set.remove((original_other_path,other_cov,other_path_type,other_path_seeds))
                            score_dict.pop((cur_path,cur_cov))
                            score_dict.pop((original_other_path,other_cov))
                            merged_in_this_iteration = True
                            break
                        else:
                            logger.info(f"bad contig merge: merged_score: {merged_score} ,skip")

                    else:
                        # Supposing that a cycle is a confident contig,namely,we randamly select a enterpoint and insert the whole path
                        

                        enterpoint = sorted(intersection, key=lambda nd: get_cov_from_spades_name(nd), reverse=True)[0]
                        logger.info(f"normal merging,  enterpoint is {enterpoint}")
                        other_score = score_dict[(other_path,other_cov)]
                        
                        merged_path = merge_paths(cur_path,other_path, enterpoint)
                        merged_score = get_score_from_path(merged_path,SEQS,max_k)

                        # 修正后的合并判断条件：
                        # 1. 评分必须不低于两者中最高的评分
                        if merged_score >= max(cur_score, other_score) * 0.95:

                            # 合并
                            logger.info(f"cur: cov: {cur_cov}, score: {cur_score}")
                            logger.info(f"other: cov: {other_cov}, score: {other_score}")
                            merged_cov  = (cur_cov + other_cov)/2
                            logger.info(f"merged_path: {merged_path} \ncov: {merged_cov}, score: {merged_score}")
                            paths_set.add((merged_path,merged_cov,str(cur_path_type+other_path_type),cur_path_seeds+other_path_seeds))
                            score_dict[(merged_path,merged_cov)] = merged_score
                            #删除原有path
                            paths_set.remove((cur_path,cur_cov,cur_path_type,cur_path_seeds))
                            paths_set.remove((original_other_path,other_cov,other_path_type,other_path_seeds))
                            score_dict.pop((cur_path,cur_cov))
                            score_dict.pop((original_other_path,other_cov))
                            merged_in_this_iteration = True
                            #直接重新进行,防止在循环中修改set可能导致的bug
                            break

                        else:
                            logger.info(f"bad merge: cur_score: {cur_score}, other_score: {other_score} --> merged_score: {merged_score}, skip")
            
            if merged_in_this_iteration:
            #直接重新进行,防止在循环中修改set可能导致的bug
                break
             
    good_paths_set = set()

    for path, cov, path_type, seeds in paths_set:
        if is_good_cyc(path,valid_mate_pairs):
            good_paths_set.add((path,cov,path_type,seeds))   
    return  good_paths_set

def get_score_from_set(paths_set: set, SEQs,max_k=77):
    score_dict = {}
    for path,cov,_ ,_ in paths_set:
        score = get_score_from_path(path,SEQs,max_k=max_k)
        score_dict[(path,cov)] = score
        # rc_p = rc_path(path)
        # if all(nd in SEQs for nd in rc_p):
        #     rc_score = get_score_from_path(rc_p,SEQs)
        #     score_dict[(rc_p,cov)] = rc_score
    return score_dict

def merge_paths(path1, path2,enterpoint):
    
    # 将短序列插入到长序列中
    if len(path1) < len(path2):
        path1 , path2 = path2, path1

    # 找到路径中交集元素的位置
    index1 = path1.index(enterpoint)
    index2 = path2.index(enterpoint)

    merged_path = path1[:index1] + path2[index2:] + path2[:index2] + path1[index1:]
    return merged_path
    
def is_similar(cur_set_mean_cov,other_set_mean_cov, cutoff=0.15,):

    # 考虑小环可能经过多次，所以小环的cov会是大环的倍数，合格后的cov应该选择大环的mean_cov
    cur_set_mean_cov = float(cur_set_mean_cov)
    other_set_mean_cov = float(other_set_mean_cov)
    # 两环mean coverage相近
    similar_difference = get_relative_difference(cur_set_mean_cov, other_set_mean_cov)
    logger.info(f"similar_diff: {similar_difference}")
    # TODO 当cov是倍数关系时，最终环的cov选择长环的，否则选两者的均值
    return  (similar_difference < cutoff) 

def get_relative_difference(a:float, b:float):
    return abs(a-b)/((a+b)/2)

def get_mean_cov_from_set(Set:set):
    total_cov = 0
    for node in Set:
        total_cov+=get_cov_from_spades_name(node)
    return total_cov/len(Set)  

def rc_path(path: list):
    # Reverse complement the nodes and create a list in one step.
    rc_path = [rc_node(node) for node in reversed(path)]  # Apply reverse and rc_node in one pass
    return tuple(rc_path)

def get_score_from_path(path,SEQS,max_k=77,num_procs=1):
    # 保存原始 stdout/stderr
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    
    # 打开 /dev/null（跨平台）
    with open(os.devnull, 'w') as devnull:
        sys.stdout = devnull
        sys.stderr = devnull
        try:
            # print("predict score for path")
            c = plasclass.plasclass()
            seqs = get_seq_from_path(path,SEQS,max_k)
            prob = c.classify(seqs)
        finally:
            # 恢复原始输出
            sys.stdout = original_stdout
            sys.stderr = original_stderr
    return prob

def get_contig_path(path_file,id_dict,SEQS,G,contig_path_file,score_out_file,min_contig_path_len=4,max_k=77,num_procs=1):
    start_time = time.time()
    # 图中所有节点的set集合
    node_set = set(G.nodes())

    # 用于存储contig path中符合条件的的path, key为设置的唯一id，value为path
    contigs_path_name_dict = {}

    # 包含两个字典，单个字典的结构为 {key:路径的首节点，value: 以key为首节点的path的list }
    # 两个字典分别考虑正向和反向序列(不是互补)，用于双向dijistra算法
    contigs_path_dict = [{},{}]

    # 存储路径对应的得分
    scores_dict = {}

    # 用于节点反向索引路径
    node_to_contig = {}

    with open(contig_path_file,'w') as o, open(path_file) as f:
        # used to identify a specific contig path with same contig name
        count = 0
        for line in f:
            if line.startswith("NODE"):
                contig_name = line.strip()
                count = 0
            else:
                line = re.sub('[:;,]'," ",line)
                path_list = line.strip().split()
                if len(path_list) >= min_contig_path_len:
                    # logger.info(f"before path_list: {str(path_list)}")
                    path_list = remove_tail_self_loop(path_list)
                    #example:  
                    #   contigs.path: 20+, 10-
                    #       ==>
                    #   assembly.fastg: EDGE_20_length_100_cov_11, EDGE_10_length_100_cov_10'
                    path_list = [transfer_to_fullname(node, id_dict) for node in path_list]

                    # trim contig path to prevent dead end
                    start_pos, end_pos = 0, len(path_list)
                    trimmed = False
                    for node in path_list:
                        # dead end
                        if node not in G.nodes():
                            trimmed = True
                            start_pos+=1
                        else:
                            break
                    for node in reversed(path_list):
                        if node not in G.nodes():
                            trimmed = True
                            end_pos-=1
                        else:
                            break
                    if trimmed:
                        logger.info(f"trimmed contig: {contig_name}:\nbefore:{path_list}\nafter:{path_list[start_pos:end_pos]} ")
                    path_list = path_list[start_pos:end_pos]

                    # Trim contig to prevent the beginning and end nodes from being the same, which may result in the inability to find a suitable path in the shortest path
                    while(len(path_list) and path_list[0] == path_list[-1]):
                        path_list = path_list[:-1]

                    if len(path_list) < min_contig_path_len or get_total_len_from_path(path_list,max_k) < 10000:
                        continue

                    # 自定义的contig path的唯一标识
                    unq_contig_name = contig_name+"_"+str(count)
                        

                    # add reversed path for bidireactional dijistra
                    freq_vec = contig_to_freq_vector(get_seq_from_path(path_list, SEQS, max_k, cycle=False))
                    contigs_path_name_dict[unq_contig_name] = (path_list, freq_vec)
                    reversed_path_list = list(reversed(path_list))
                    contigs_path_name_dict['R'+unq_contig_name] = (reversed_path_list, freq_vec)
                    count+=1
                    # store contig path 
                    # the judgement may be redundant
                    if all(node in node_set for node in path_list ):
                        contigs_path_dict[0].setdefault(path_list[0],[])
                        contigs_path_dict[1].setdefault(reversed_path_list[0],[])

                        # dont store the first node which stored in key
                        # key:  first node of contig path
                        # value: (contig path, contig path name, score of contig path, conitg path fragment end with key(here is Null))
                        contigs_path_dict[0][path_list[0]].append((path_list[1:],unq_contig_name,0,None,freq_vec))
                        contigs_path_dict[1][reversed_path_list[0]].append((reversed_path_list[1:],'R'+unq_contig_name,0,None,freq_vec))
                            
        logger.info(f"preprocess contig path score...")
        print("Preprocessing contig path score")
        

        # 将所有的正向contig path写入文件并批量预测得分
        for node, info in contigs_path_dict[0].items():
            for path, name,_,_,freq_vec in info:
                if len(path) >= min_contig_path_len-1:
                    o.write(f">{name}\n")
                    new_path = [node]+path[:]
                    o.write(f"{get_seq_from_path(new_path,SEQS,max_k,False)}\n")  
                else:
                    raise ValueError(f"path of {name} : {path} have illeagel length")
    # 保存原始 stdout/stderr
    original_stdout = sys.stdout
    original_stderr = sys.stderr

    # 打开 /dev/null（跨平台）
    with open(os.devnull, 'w') as devnull:
        sys.stdout = devnull
        sys.stderr = devnull
        try:
            classify(contig_path_file, score_out_file, num_procs)
        finally:
            # 恢复原始输出
            sys.stdout = original_stdout
            sys.stderr = original_stderr
    # 存储得分, plasclass 的输出格式： contig_name      score
    with open(score_out_file,'r')as f:
        for line in f:
            cname, score = line.strip().split()[:2]
            score  = transformByLength([get_total_len_from_path(path,max_k)], [float(score)])[0]
            scores_dict[cname] = score
           
    
    #将得分更新入字典中，并筛除得分较低的contig path
    for k in (0, 1):
        for node, info_list in list(contigs_path_dict[k].items()):  # 使用list()创建副本以安全地修改原字典
            to_remove = []  # 创建一个列表来存储需要移除的索引
            for idx, (path, name, s, pre,freq_vec) in enumerate(info_list):
                if len(path) >= min_contig_path_len-1:
                    # 因为仅计算了正向序列的得分，所以反向序列需要去掉前缀'R'来获取正确的名称
                    score = scores_dict.get(name, None) if k == 0 else scores_dict.get(name[1:], None)
                    if score is not None and score >= 0.5:
                        info_list[idx] = (path, name, score, None,freq_vec)  # 更新整个元组
                        for n in path:      #更新节点到path的索引
                            node_to_contig.setdefault(n,set())
                            node_to_contig[n].add(name)
                    else:
                        to_remove.append(idx)  # 记录要删除的元组索引
                else:
                    raise ValueError(f"path of {node} has illegal length: {len(path)}")
            
            # 删除不符合条件的元组，逆序删除以避免索引问题
            for idx in reversed(to_remove):
                del info_list[idx]

        end_time = time.time()
    logger.info(f"preprocess finished, consuming {end_time-start_time} seconds")

    return contigs_path_dict ,contigs_path_name_dict, node_to_contig, scores_dict

def transfer_to_fullname(node:str, id_dict):
    if node.endswith('+'):
        return id_dict[node[:-1]]
    elif node.endswith('-'):
        return id_dict[node[:-1]]+"'"
    else:
        sys.exit(f"contig path file format error: {node} must end with '+' or '-' ")

def dijkstra_path(G, path_dict:dict,node_score_dict,node_vec_dict,node_support_dict,SEQS, source, target, weight='weight',bidirectional=False):


    (length, path,use_contig) = single_source_dijkstra(G, path_dict,node_score_dict,node_vec_dict,node_support_dict,SEQS ,source, target=target,weight=weight,bidirectional=bidirectional)
    return length,path,use_contig

def single_source_dijkstra(G, path_dict,node_score_dict,node_vec_dict,node_support_dict,SEQS,source, target=None, cutoff=None,
                           weight='weight',bidirectional=False,max_k_val=77):
    use_contig = None
    if not source:
        raise ValueError('sources must not be empty')
    if target == source:
        return (0, [target],[])

    weight_lambda = _weight_function(G, weight)
    G_succ = G._succ if G.is_directed() else G._adj
    node_set =set(G_succ.keys())
    min_len = float("inf")
    final_path = None
    target_info = (None,None,None)
    if target is None:
        return (min_len, [source],use_contig)
   
    # 单向dijstra只利用path_dict中存储的正向的contig path,即path_dict[0]
    if source in path_dict[0].keys():

        for record in path_dict[0][source]:
            try:
                path,contig_name,score,pre_contig,freq_vec = record
                
                # logger.info(f"path dict has {contig_name}")
                finaldist, finalpath = float("inf"),None
                # 当前节点为contig path的首节点,直接找contig path的尾节点到当前节点的前驱节点的最短路径即可
                if pre_contig is None:
                    cur_path = {path[-1] : [source]+path}
                    if all(node in node_set for node in path):
                        
                        finaldist, finalpath, target_info = _dijkstra_multisource(G, path_dict,node_score_dict,node_vec_dict,node_support_dict,SEQS,[[source]+path,[score, freq_vec] ], weight_lambda, paths=cur_path,
                                    cutoff=cutoff, target=target)
                        # else:
                        #     finaldist, finalpath,_, target_info = bidirectional_dijkstra(G, path_dict,SEQS,source, target,weight,max_k_val=max_k_val)
                        target, target_record, used_contigs = target_info
                        if target_record is not None:
                            u_score, u_freq_vec = target_record
                        else:
                            assert len(target) == 1
                            u_score, u_freq_vec = node_score_dict[canonicalize(target[0])],node_vec_dict[canonicalize(target[0])]
                        finaldist+= get_edge_cost(G, target, [source]+path, u_score, score, u_freq_vec, freq_vec,node_support_dict.get((target[-1],source),0), max_k_val)
                        # add first contig to final path

                    else:
                        # logger.info(f"{contig_name} not in G: {', '.join(str(node) for node in path if node not in node_set)}")
                        continue

                if finaldist < min_len:
                    min_len = finaldist
                    final_path = finalpath
                    # logger.info(f"final add contig path: {sources[0]}--->{target}: {finalpath}, weight: {finaldist}")
                    use_contig = [contig_name] + target_info[2]
            except KeyError:
                # logger.info(f"has not path from {sources[0]} to {target} through {contig_name}")
                continue
    # can't find path through contig
    # 如果通过contig找不到最短环,则直接通过双向dijistra找
    if final_path is None:
        # finaldist, finalpath, through_contig= bidirectional_dijkstra(G, path_dict,SEQS,source, target,weight,max_k_val=max_k_val)
        cur_path = {source:[source]}
        finaldist, finalpath, target_info= _dijkstra_multisource(G, path_dict,node_score_dict,node_vec_dict,node_support_dict,SEQS,[[source],None], weight_lambda, paths=cur_path,
                                 cutoff=cutoff, target=target)
        target, target_record, _ = target_info
        if target_record is not None:
            u_score, u_freq_vec = target_record
        else:
            assert len(target) == 1
            u_score, u_freq_vec = node_score_dict[canonicalize(target[0])], node_vec_dict[canonicalize(target[0])]
        finaldist+= get_edge_cost(G, target, [source], u_score, node_score_dict[canonicalize(source)], u_freq_vec, node_vec_dict[canonicalize(source)],node_support_dict.get((target[-1],source),0), max_k_val)
        try:
            # logger.info(f"simple path: {sources[0]}--->{target}: {finalpath}, weight: {finaldist}")
            min_len = finaldist
            final_path = finalpath
            use_contig = target_info[2]
        except KeyError:
            raise nx.NetworkXNoPath("No path to {}.".format(target))

    return (min_len,final_path,use_contig)  
def _weight_function(G, weight):

    if callable(weight):
        return weight
    # If the weight keyword argument is not callable, we assume it is a
    # string representing the edge attribute containing the weight of
    # the edge.
    if G.is_multigraph():
        return lambda u, v, d: min(attr.get(weight, 1) for attr in d.values())
    return lambda u, v, data: data.get(weight, 1)


def _dijkstra_multisource(G, path_dict: dict,node_score_dict,node_vec_dict,node_support_dict, SEQS, source, weight, pred=None, paths=None,
                                cutoff=None, target=None, max_k_val=77):
    """
    使用 Min-Max 松弛规则的 Dijkstra 风格路径搜索，寻找瓶颈最小的路径。
    在瓶颈权重相同时，使用路径长度（边数）作为平局打破规则，选择更短的路径。
    """
    G_succ = G._succ if G.is_directed() else G._adj
    node_set = set(G_succ.keys())
    push = heappush
    pop = heappop
    
    # dist 存储 Min-Max 权重（瓶颈 Cost）
    dist = {}
    # seen 存储在 fringe 中的 Min-Max 权重
    seen = {}
    c = count()
    fringe = []
    target_info = (None, None, None)

    # source 应为 (node_list, record)
    start_nodes = source[0]
    start_record = source[1] if len(source) > 1 else None
    
    start_node = start_nodes[-1]
    if start_node not in G:
        raise nx.NodeNotFound(f"Source {start_node} not in G")
    
    # 路径长度 (Length): 初始路径长度为 0
    start_length = len(start_nodes) - 1 if start_nodes else 0
    
    # Min-Max Cost: 源点到自身的瓶颈 Cost 为 0
    start_cost = 0 
    seen[start_node] = (start_cost, start_length)
    
    # 优先队列元素格式: (MinMax_Cost, Length, Tie_Breaker, (Path, Record, Used_Contigs))
    push(fringe, (start_cost, start_length, next(c), (start_nodes, start_record, [])))

    while fringe:
        (d, length, _, (v, cur_record, used_contigs)) = pop(fringe)
        
        # v[-1] 是当前路径的终点节点
        if v[-1] in dist:
            continue
            
        # 确认 Min-Max 路径
        dist[v[-1]] = d

        if v[-1] == target:
            # 找到目标，返回信息
            target_info = (v, cur_record, used_contigs)
            break

        # ----------------------------------------------------
        # 1. 松弛操作：尝试走 contig path (Contig-Path-Contig)
        # ----------------------------------------------------
        through_contig = False
        for u_start_node, e in G_succ[v[-1]].items():
            if u_start_node in path_dict[0]:
                for record in path_dict[0][u_start_node]:
                    u_path, u_contig_name, u_score, u_pre_contig, u_vec = record
                    
                    
                    # 获取节点评分/向量信息
                    v_score, v_vec = cur_record if cur_record is not None else (node_score_dict[canonicalize(v[0])], node_vec_dict[canonicalize(v[0])])
                    
                    # 边 (v, Contig Path) 的 Cost (瓶颈权重)
                    cost = get_edge_cost(G, v, [u_start_node] + u_path, v_score, u_score, v_vec, u_vec,node_support_dict.get((v[-1], u_start_node),0), max_k_val)
                    
                    # **Min-Max 松弛规则**
                    # vu_dist = max(dist[v[-1]], cost) 
                    vu_dist = dist[v[-1]] + cost
                    
                    # 新路径长度 (每经过一个 Contig 视为增加一个边)
                    new_length = length + len(u_path) 
                    
                    u_end_node = u_path[-1] # Contig Path 的终点

                    # 获取目标节点当前的 Min-Max 信息
                    current_minmax_cost, current_length = seen.get(u_end_node, (float('inf'), float('inf')))

                    # 比较规则: 
                    # 1. 新的 Min-Max 权重必须更小 (vu_dist < current_minmax_cost)
                    # 2. 如果 Min-Max 权重相同，新的路径必须更短 (vu_dist == current_minmax_cost AND new_length < current_length)
                    if vu_dist < current_minmax_cost or \
                        (vu_dist == current_minmax_cost and new_length < current_length):
                        through_contig = True
                        # 更新 seen 和 priority queue
                        seen[u_end_node] = (vu_dist, new_length)
                        new_contigs = used_contigs + [u_contig_name]
                        
                        # 将 length 作为第二个排序元素 (最小化)
                        push(fringe, (vu_dist, new_length, next(c), (u_path, [u_score, u_vec], new_contigs)))
                        
                        if paths is not None:
                            paths[u_end_node] = paths.get(v[-1], v) + [u_start_node] + u_path


        # ----------------------------------------------------
        # 2. 松弛操作：走普通边 (Contig-Contig)
        # ----------------------------------------------------
        if not through_contig:
            for u, e in G_succ[v[-1]].items():
                u_score, u_vec = node_score_dict[canonicalize(u)], node_vec_dict[canonicalize(u)]
                v_score, v_vec = cur_record if cur_record is not None else (node_score_dict[canonicalize(v[0])], node_vec_dict[canonicalize(v[0])])
                
                # 边 (v, u) 的 Cost (瓶颈权重)
                cost = get_edge_cost_from_graph(G, v, [u], v_score, u_score, v_vec, u_vec,node_support_dict.get((v[-1],u),0), max_k_val)
                
                # **Min-Max 松弛规则**
                # vu_dist = max(dist[v[-1]], cost) 
                vu_dist = dist[v[-1]] + cost
                
                # 新路径长度
                new_length = length + 1

                # 获取目标节点当前的 Min-Max 信息
                current_minmax_cost, current_length = seen.get(u, (float('inf'), float('inf')))
                
                # 比较规则 (同 Contig Path 比较)
                if vu_dist < current_minmax_cost or \
                   (vu_dist == current_minmax_cost and new_length < current_length):
                    
                    # 更新 seen 和 priority queue
                    seen[u] = (vu_dist, new_length)
                    push(fringe, (vu_dist, new_length, next(c), ([u], None, used_contigs)))
                    
                    if paths is not None:
                        paths[u] = paths.get(v[-1], v) + [u]

    else:
        # 如果循环结束仍未找到目标
        raise nx.NetworkXNoPath(f"No path between {source[0]} and {target}.")

    finaldist = dist[target]
    finalpath = paths[target] if paths else None
    return finaldist, finalpath, target_info

def classify(infile, outfile, num_procs):
    ''' Run the classification
    '''
    c = plasclass.plasclass(num_procs)
    seq_names = []
    seqs = []
    i = 0
    fp = open(infile)
    with open(outfile,'w') as o:
        for name, seq, _ in readfq(fp):
            seq_names.append(name)
            seqs.append(seq)
            i += 1
            if i % 50000 == 0:
                probs = c.classify(seqs)
                for j,p in enumerate(probs):
                    o.write(seq_names[j] + '\t' + str(p) + '\n')
                seq_names = []
                seqs = []


        # last bunch of sequences:
        probs = c.classify(seqs)
        for j,p in enumerate(probs):
            o.write(seq_names[j] + '\t' + str(p) + '\n')

    fp.close()

def get_total_len_from_path(path, max_k_val, cycle=False):
    total_len = 0
    total_len = sum(get_length_from_spades_name(nd) - max_k_val for nd in path)
    if cycle is False:
        total_len += max_k_val
    return total_len


def add_contig_to_path_dict(
    G, scores_dict, node_score_dict,node_gene_set,path_dict, contigs_path_name_dict,
    node_to_contig, use_genes=True, use_scores=True
):
    # Step 1: 预计算每个 contig 中节点的位置
    proxy_contig_dict = {}
    contig_node_to_index = {}
    for contig_id, (path, _) in contigs_path_name_dict.items():
        contig_node_to_index[contig_id] = {node: i for i, node in enumerate(path)}

    # Step 2: 确定高置信节点
    plasmid_nodes = set()
    if use_genes:
        plasmid_nodes.update(get_plasmid_gene_nodes(G,node_gene_set))
    if use_scores:
        plasmid_nodes.update(get_hi_conf_plasmids(G,node_score_dict))

    # Step 3: 遍历每个高置信节点
    for node in plasmid_nodes:
        if node not in node_to_contig:
            continue

        for contig_id in node_to_contig[node]:
            if contig_id.startswith('R'):  # skip reverse
                continue

            # 快速获取该 node 在 contig 中的位置
            node_index_map = contig_node_to_index.get(contig_id)
            if node_index_map is None or node not in node_index_map:
                continue

            start_index = node_index_map[node]
            path, freq_vec = contigs_path_name_dict[contig_id]

            # Skip if at boundary
            if start_index == 0 or start_index == len(path) - 1:
                continue

            # Split path
            contig_path1 = path[start_index:]      # [node, ..., end]
            contig_path2 = path[:start_index]      # [start, ..., prev]

            score = scores_dict.get(contig_id)
            if score is None:
                continue

            # Store: key = first node of suffix path
            suffix_start = contig_path1[0]
            suffix_body = contig_path1[1:]  # exclude first node

            proxy_contig_dict.setdefault(suffix_start, []).append(
                (suffix_body, contig_id, score, contig_path2, freq_vec)
            )

            #     logger.info(f"hi conf node: {node}")
            #     logger.info(f"add contig path: start={suffix_start}, body={suffix_body}, "
            #                 f"score={score}, id={contig_id}, pre={contig_path2}")    
    return proxy_contig_dict 
        
             
def remove_dead_ends(G):

    """
    Recursively removes all nodes with zero in-degree or out-degree from a directed graph.

    Parameters:
    G (networkx.DiGraph): The directed graph to process.

    """
    has_dead_end = True

    while has_dead_end:
        # Collect nodes to remove in this iteration
        nodes_to_remove = [node for node in G if G.in_degree(node) == 0 or G.out_degree(node) == 0]
        
        if not nodes_to_remove:
            has_dead_end = False
        else:
            # Remove all collected nodes at once
            G.remove_nodes_from(nodes_to_remove)

def merge_path_through_contig_path(cur_path, other_path, nd, node_to_contig, contigs_path_name_dict):
    # 获取 nd 所属的所有 contig ID
    # logger.info(f"tempt to merge:\n{str(cur_path)}\n{str(other_path)}\n through {nd}\n")
    
    contig_ids = node_to_contig.get(nd, [])
    
    for contig_id in contig_ids:
        contig = contigs_path_name_dict.get(contig_id, None)
        if contig is None or len(contig) < 3:
            continue  # 忽略不存在或长度不足的 contig
        
        try:
            idx = contig.index(nd)
        except ValueError:
            continue  # 忽略 nd 不在 contig 中的情况

        # 处理边界情况
        if idx == 0 or idx == len(contig) - 1:
            continue  # 忽略 nd 是 contig 的第一个或最后一个元素的情况
        # logger.info(f"through contig {str(contig)}")
        min_shared_part = (contig[idx-1], nd, contig[idx+1])
        # logger.info(f"shared part: {str(min_shared_part)}")
        # 串联形成头尾结合的部分
        double_cur_path = cur_path+cur_path
        double_other_path = other_path+other_path

        left_part_in_cur_path = contain(double_cur_path,min_shared_part[:2])
        right_part_in_other_path = contain(double_other_path,min_shared_part[1:])
        # logger.info(f"left_part_in_cur_path:{left_part_in_cur_path}, right_part_in_other_path: {right_part_in_other_path}")

        if left_part_in_cur_path != -1 and right_part_in_other_path!= -1:
            return merge_contig_path(cur_path, other_path, nd, left_part_in_cur_path,right_part_in_other_path)
        

        left_part_in_other_path = contain(double_other_path,min_shared_part[:2])
        right_part_in_cur_path = contain(double_cur_path,min_shared_part[1:])
        # logger.info(f"left_part_in_other_path:{left_part_in_other_path}, right_part_in_cur_path: {right_part_in_cur_path}")

        if left_part_in_other_path != -1 and right_part_in_cur_path!= -1:
            return merge_contig_path(other_path, cur_path, nd, left_part_in_other_path,right_part_in_cur_path)
           
    
    return None

def contain(path, subseq):
    # logger.info(f"check {subseq} in \n {path}")
    for i in range(len(path)//2):
        if path[i:i+len(subseq)] == subseq:
            return i+1
    return -1

def merge_contig_path(path1, path2, nd,left_idx, right_idx):
    # logger.info(f"merge {str(path1)} and {str(path2)}")
    right_idx-=1
    return path1[:left_idx] + path2[right_idx:]+path2[:right_idx] + path1[left_idx:]

def meet_criterion(path, G, SEQS,max_k, max_CV,valid_pairs):
    return get_wgtd_path_coverage_CV(path,G,SEQS,max_k_val=max_k) <= max_CV and is_good_cyc(path,valid_pairs)

def sort_key(path):
    return not path[2]

def extract_node_id(node_label):
    """
    从节点标签中提取出ID部分,格式为 "EDGE_<id>_length_<length>_cov_<coverage>"
    """
    reverse = False
    mark = None
    match = re.match(r"EDGE_(\d+)_length_\d+_cov_[\d.]+", node_label)
    if not match:
        raise ValueError(f"Unexpected node label format: {node_label}")
    if node_label[-1] == "'":
        reverse = True
    if reverse:
        mark = match.group(1)+"-"
    else:
        mark = match.group(1)+"+"
    return mark

def component_sort_key(comp):
    """
    创建一个排序键, 基于组件中所有节点的ID拼接而成的字符串。
    """
    # 提取并排序节点ID
    sorted_node_ids = sorted((extract_node_id(n) for n in comp.nodes()), key=str)
    
    # 将排序后的节点ID拼接成字符串
    id_string = ','.join(str(id) for id in sorted_node_ids)
    
    return id_string

def remove_tail_self_loop(path):
    """
    剔除路径尾部与开头重复的部分（模拟自环冗余）
    
    输入: path = ['71879+', '72249+', '71881+', '72249+', '72067-', '71643-', '71879+', '72249+']
    输出: ['71879+', '72249+', '71881+', '72249+', '72067-', '71643-']
    """
    if not path or len(path) < 2:
        return path

    # Step 1: 构建 KMP 的 LPS 数组（最长公共前后缀长度）
    def compute_lps(pattern):
        lps = [0] * len(pattern)
        length = 0  # 当前最长公共前后缀的长度
        i = 1
        while i < len(pattern):
            if pattern[i] == pattern[length]:
                length += 1
                lps[i] = length
                i += 1
            else:
                if length != 0:
                    length = lps[length - 1]
                else:
                    lps[i] = 0
                    i += 1
        return lps

    lps = compute_lps(path)
    
    # lps[-1] 是整个序列的最长公共前后缀长度
    overlap_len = lps[-1]
    
    # 如果 overlap_len > 0，说明末尾 overlap_len 个元素等于开头 overlap_len 个元素
    # 我们将其从尾部删除（保留开头，去掉尾部重复）
    if overlap_len > 0:
        # 检查是否真的是尾部匹配开头
        if path[-overlap_len:] == path[:overlap_len]:
            # 删除尾部 overlap_len 个元素
            return path[:-overlap_len]
    
    return path

def get_edge_cost_from_graph(G, from_path, to_path, from_score, to_score, from_vec, to_vec,support, max_k_val):
    # normal edge
    if(len(from_path) == 1 and len(to_path) == 1):
        return G[from_path[0]][to_path[0]]['cost']
    # through contig path
    else:
        return get_edge_cost(G, from_path, to_path, from_score, to_score, from_vec, to_vec,support, max_k_val)
    
def get_native_path_dict(COMP, contigs_path_dict):
    comp_nodes = set(COMP.nodes())  
    path_dict = [{}, {}]

    for node, info_list in contigs_path_dict[0].items():
        if node not in comp_nodes:
            continue
        filtered = []
        for path, name, s, pre, freq_vec in info_list:
            if all(nd in comp_nodes for nd in path):  
                filtered.append((path, name, s, pre, freq_vec))
        if filtered:
            path_dict[0][node] = filtered
    return path_dict


def get_native_proxy_path_dict(COMP, all_proxy_path_dict):
    comp_nodes = set(COMP.nodes()) 
    proxy_path_dict = {}
    for node, info_list in all_proxy_path_dict.items():
        if node not in comp_nodes:
            continue
        filtered = []
        for path, name, s, pre, freq_vec in info_list:
            if all(nd in comp_nodes for nd in path):
                filtered.append((path, name, s, pre, freq_vec))
        if filtered:
            proxy_path_dict[node] = filtered
    return proxy_path_dict

def trim_path_tips_light(G, path):
    if len(path) <= 1:
        return path
    
    # 只关注 path 内部的连接（构建 induced subgraph 的度数）
    path_set = set(path)
    in_deg = {}
    out_deg = {}
    
    for v in path:
        # 只统计 path 内部的邻居
        in_deg[v] = sum(1 for u in G.predecessors(v) if u in path_set)
        out_deg[v] = sum(1 for w in G.successors(v) if w in path_set)
    
    path = list(path)
    i, j = 0, len(path) - 1
    
    # 从前向后 trim sources
    while i < j and in_deg[path[i]] == 0:
        u = path[i]
        # 更新后继的入度
        for w in G.successors(u):
            if w in path_set and i < path.index(w) <= j:
                in_deg[w] -= 1
        i += 1
    
    # 从后向前 trim sinks
    while i < j and out_deg[path[j]] == 0:
        u = path[j]
        for w in G.predecessors(u):
            if w in path_set and i <= path.index(w) < j:
                out_deg[w] -= 1
        j -= 1
    
    return path[i:j+1]

def canonicalize(nd: str):
    return nd[:-1] if nd.endswith("'") else nd

def get_shortest_batch(batch_nodes):

    G = worker_data['G']
    SEQS = worker_data['SEQS']
    node_score_dict = worker_data['scores']
    path_dict = worker_data['path']
    node_vec_dict = worker_data['vec']
    node_support_dict = worker_data['support']

    all_paths = []
    for node in batch_nodes:
        res_list = get_shortest((node, path_dict, node_score_dict, node_vec_dict, node_support_dict, SEQS, G))
        if res_list:
            for path_list, weight, use_contig in res_list:
                all_paths.append((tuple(path_list), weight, use_contig,node))  
    return all_paths


def calc_anchored_cov(node, path, G, max_depth=20):
    """
    计算节点的折算覆盖度 (最小流量法)。
    核心思想：路径的覆盖度受限于上下游最窄的瓶颈（木桶效应）。
    """
    if node not in G: return 0.0
    
    node_cov = get_cov_from_spades_name_and_graph(node, G) # 适配你的 get_cov 函数
    path_set = set(path) # 假设已处理 RC
    
    # --- 辅助函数：计算某节点在特定方向上的有效路径流量 ---
    def get_path_flow(anchor, direction):
        if direction == 'pred':
            # 作为下游，看上游(Predecessors)有多少流量进来
            neighbors = list(G.predecessors(anchor))
        else:
            # 作为上游，看下游(Successors)有多少流量流出
            neighbors = list(G.successors(anchor))
            
        total_cov = sum(get_cov_from_spades_name_and_graph(n, G) for n in neighbors)
        if total_cov <= 1e-9: return 0.0
        
        # 只累加流向/来自 Path 的邻居的覆盖度
        path_flow = sum(get_cov_from_spades_name_and_graph(n, G) for n in neighbors if n in path_set)
        
        # 这里有一个关键点：
        # 如果是单纯的线性连接，Flow = Neighbor_Cov。
        # 如果是分叉/汇合，Flow = 分支流量。
        return path_flow

    # ==========================================
    # 1. 向上游寻找流量瓶颈 (Upstream Flow)
    # ==========================================
    upstream_flow = float('inf')
    curr = node
    
    # 如果当前节点就是 Contig 起点(无前驱)，那上游流量就是自身
    if G.in_degree(curr) == 0:
        upstream_flow = node_cov
    else:
        for _ in range(max_depth):
            preds = list(G.predecessors(curr))
            
            # 停止条件 A: 遇到汇合点 (Merge Point)
            # 此时我们要看这个汇合点给当前路径贡献了多少
            if len(preds) != 1:
                upstream_flow = get_path_flow(curr, 'pred')
                break
                
            prev = preds[0]
            
            # 停止条件 B: 上一个是分叉点 (Split Point)
            # 意味着上游在分流，我们要看分给我们的这一支有多少
            if G.out_degree(prev) != 1:
                # 注意：这里我们计算 prev 的 output flow (succs)，
                # 理论上应该等于 curr 的 input flow。
                # 直接用 curr 的 input flow 计算更简单且等价。
                upstream_flow = get_path_flow(curr, 'pred')
                break
            
            # 停止条件 C: 路径断裂
            if prev not in path_set:
                upstream_flow = get_path_flow(curr, 'pred')
                break
                
            curr = prev
        else:
            # 超过深度，默认信任当前节点的输入
            upstream_flow = get_path_flow(node, 'pred')

    # ==========================================
    # 2. 向下游寻找流量瓶颈 (Downstream Flow)
    # ==========================================
    downstream_flow = float('inf')
    curr = node
    
    if G.out_degree(curr) == 0:
        downstream_flow = node_cov
    else:
        for _ in range(max_depth):
            succs = list(G.successors(curr))
            
            # 停止条件 A: 遇到分叉点 (Split Point)
            # 我们要看这个点分给路径多少
            if len(succs) != 1:
                downstream_flow = get_path_flow(curr, 'succ')
                break
            
            next_n = succs[0]
            
            # 停止条件 B: 下一个是汇合点 (Merge Point)
            # 我们流向了一个大熔炉，我们的贡献由当前的输出决定
            if G.in_degree(next_n) != 1:
                downstream_flow = get_path_flow(curr, 'succ')
                break
                
            if next_n not in path_set:
                downstream_flow = get_path_flow(curr, 'succ')
                break
                
            curr = next_n
        else:
            downstream_flow = get_path_flow(node, 'succ')

    # ==========================================
    # 3. 取最小值 (Min-Flow)
    # ==========================================
    # 逻辑：一条管子的最大流量取决于它最窄的地方
    
    # 修正 inf (如果未找到有效上/下游，说明没有瓶颈，取自身)
    if upstream_flow == float('inf'): upstream_flow = node_cov
    if downstream_flow == float('inf'): downstream_flow = node_cov
    
    return min(node_cov, upstream_flow, downstream_flow)
