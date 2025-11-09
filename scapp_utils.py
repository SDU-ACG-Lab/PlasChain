from __future__ import division
import numpy as np
import math,copy
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

import PARAMS

complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
logger = logging.getLogger("scapp_logger")
path_find_time_consume = 0

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


def get_node_scores(scores_file,G):
    """ Write the plasmid scores into each node in the graph
    """
    scores = {}
    with open(scores_file) as f:
        for line in f:
            split = line.strip().split()
            scores[split[0]] = float(split[1])
    for nd in G.nodes():
        G.add_node(nd, score=scores[nd])

def get_gene_nodes(genes_file,G):
    """ Annotate each node in the graph whether it has plasmid gene in it
    """
    gene_nodes = set()
    with open(genes_file) as f:
        for line in f:
            gene_nodes.add(line.strip())
    for nd in G.nodes():
        if nd in gene_nodes:
            G.add_node(nd, gene=True)
        else:
            G.add_node(nd,gene=False)

def rc_seq(dna):
    rev = reversed(dna)
    return "".join([complements[i] for i in rev])

# def get_num_from_spades_name(name):
#     name_parts = name.split("_")
#     contig_length = name_parts[1]
#     return int(contig_length)
    
def get_num_from_spades_name(name):
    name_parts = name.split("_")
    contig_length = name_parts[1]
    return contig_length

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
    if new_cov == 0:
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
    if coverage <= 0.0: coverage = 1.0/float(length) # ensure no division by zero, consider more principled way to do this
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


def get_path_covs(path,G,discount=False):
    cnts = {}
    if discount:
        covs = [get_discounted_node_cov(n,path,G) for n in path]
        # discount weight of nodes that path passes through multiple times
        cnts = get_node_cnts_hist(path)
        for i in range(len(path)):
            p = path[i]
            pos_name = p if (p[-1]!="'") else p[:-1]
            if cnts[pos_name] > 1:
                covs[i] /= cnts[pos_name]
    else:
        covs = [get_cov_from_spades_name_and_graph(n,G) for n in path]

    return covs

def get_path_mean_std(path, G, seqs, max_k_val=77,discount=True):
    covs = np.array(get_path_covs(path,G,discount))
    wgts = np.array([(get_length_from_spades_name(n)-max_k_val) for n in path])
    tot_len = len(get_seq_from_path(path, seqs, max_k_val, cycle=True))
    if tot_len<=0: return (0,0)
    wgts = np.multiply(wgts, 1./tot_len)
    mean = np.average(covs, weights = wgts)
    std = np.sqrt(np.dot(wgts,(covs-mean)**2))
    return (mean,std)

def update_path_coverage_vals(path, G, seqs, max_k_val=77,):
    mean, _ = get_path_mean_std(path, G, seqs, max_k_val) ## NOTE: CAN WE STILL GUARANTEE CONVERGENCE WHEN DISCOUNTING COVERAGE ??!
    covs = get_path_covs(path,G)
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
    logger.info("Path: %s \nMean: %s Covs: %s" % (str(path),str(mean),str(covs) ))
    logger.info("NewCovs: %s" % (str(new_covs)))
    for i in range(len(path)):
        if new_covs[i] > 0:
            update_node_coverage(G,path[i],new_covs[i])
        else:
            update_node_coverage(G,path[i],0)
    return new_covs

def update_path_with_covs(path, G, covs):
    for i in range(len(path)):
        if covs[i] > 0:
            update_node_coverage(G,path[i],covs[i])
        else:
            update_node_coverage(G,path[i],0)

def get_total_path_mass(path,G):
    return sum([get_length_from_spades_name(p) * \
        get_cov_from_spades_name_and_graph(p,G) for p in path])


def get_long_self_loops(G, min_length, seqs, bamfile, use_scores=True, use_genes=True, max_k_val=77, score_thresh=0.9, mate_thresh = 0.1):
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
                if G.nodes[nd]['score'] > PARAMS.SELF_LOOP_SCORE_THRESH or G.nodes[nd]['gene']==True:
                    potential_plasmids.add(nd_path)
                    logger.info("Added path: %s - high scoring long self-loop" % nd)
                    to_remove.append(nd)
                    continue

            off_node_mate_count, on_node_mate_count = count_selfloop_mates(nd,bamfile,G)
            if float(off_node_mate_count) > PARAMS.SELF_LOOP_MATE_THRESH*float(on_node_mate_count):
                logger.info('Self loop %s has %2f percent off-node mate-pairs. Removing' % (nd,PARAMS.SELF_LOOP_MATE_THRESH))
                to_remove.append(nd)
            else:
                potential_plasmids.add(nd_path)
                logger.info("Added path: %s  - long self loop" % nd)
                to_remove.append(nd)
        else: # non-isolated loop
            if path_len < min_length: continue

            off_node_mate_count, on_node_mate_count = count_selfloop_mates(nd,bamfile,G)
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

def remove_hi_confidence_chromosome(G):
    """ Remove the long nodes that are predicted to likely be chromosomal
    """
    # dont remove chromosome nodes
    # return
    to_remove = []
    for nd in G.nodes():
        if get_length_from_spades_name(nd) > PARAMS.CHROMOSOME_LEN_THRESH and \
            G.nodes[nd]['score'] < PARAMS.CHROMOSOME_SCORE_THRESH:
            to_remove.append(nd)

            to_remove.append(rc_node(nd))
    G.remove_nodes_from(to_remove)
    logger.info(f"to remove: {str(to_remove)}")
    logger.info("Removed %d long, likely chromosomal nodes" % len(set(to_remove)))

def get_hi_conf_plasmids(G):
    """ Return a list of nodes that are likely plasmids
    """

    hi_conf_plasmids = [nd for nd in G.nodes() if (get_length_from_spades_name(nd) > PARAMS.PLASMID_LEN_THRESH and \
                        G.nodes[nd]['score'] > PARAMS.PLASMID_SCORE_THRESH)]

    logger.info("Found %d long, likely plasmid nodes" % len(hi_conf_plasmids))
    return hi_conf_plasmids

def get_plasmid_gene_nodes(G):
    """ Return list of nodes annotated as having a plasmid gene
    """
    plasmid_gene_nodes = [nd for nd in G.nodes() if G.nodes[nd]['gene']==True]
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

def get_shortest(args_array):
    """ Worker function for getting shortest path to each node in parallel
    """

    node, path_dict,SEQS,G, paths_list = args_array
    shortest_score = float("inf")
    path = None
    use_contig = False


    for pred in G.predecessors(node):
        try:
            #  path_len,shortest_path = nx.bidirectional_dijkstra(G, node, pred, weight='cost')
            # path_len,shortest_path = bidirectional_dijkstra(G,path_dict,SEQS,source=node,target=pred, weight='cost')
            path_len,shortest_path,use = dijkstra_path(G,path_dict,SEQS,source=node,target=pred, weight='cost',bidirectional=True)
            
            # logger.info(f"shortest_path:  {str(folded_path)}")
            if path_len < shortest_score:
                path = shortest_path
                shortest_score = path_len
                use_contig = use
        except nx.exception.NetworkXNoPath:
            continue
    # logger.info(f"add shortest path: {path}")
    if path is not None: paths_list.append((path,shortest_score,use_contig))
    # done

def enum_high_mass_shortest_paths(G, pool,path_dict,SEQS, use_scores=False, use_genes=False, seen_paths=None):
    """ given component subgraph, returns list of paths that
        - is non-redundant (includes) no repeats of same cycle
        - includes shortest paths starting at each node n (assigning
        node weights to be 1/(length * coverage)) to each of
        its predecessors, and returning to n
    """
    if seen_paths == None:
        seen_paths = []
    unq_sorted_paths = set([])
    # in case orientation obliv. sorted path strings passed in
    for p in seen_paths:
        unq_sorted_paths.add(p)
    paths = []

    nodes = []
    nodes = list(G.nodes()) # creates a copy

    # logger.info("Getting edge weights")

    # use add_edge to assign edge weights to be 1/mass of starting node
    # TODO: only calculate these if they haven't been/need to be updated
    for e in G.edges():
        if use_genes and G.nodes[e[1]]['gene'] == True:
            G.add_edge(e[0], e[1], cost = 0.0)
        elif use_scores==True:
            G.add_edge(e[0], e[1], cost = (1.-(G.nodes[e[1]]['score']))/get_spades_base_mass(G, e[1]))
        else:
            G.add_edge(e[0], e[1], cost = (1./get_spades_base_mass(G, e[1])))

    logger.info("Getting shortest paths")
    st = time.time()
    paths_list = []
    if pool._processes > 1 and pool._processes <= 2*len(nodes): # otherwise, run single threaded
        paths_list=Manager().list()
        pool.map(get_shortest, [[node, path_dict,SEQS,G, paths_list] for node in nodes])
    else:
        for node in nodes:
            get_shortest([node,path_dict,SEQS,G,paths_list])
    ed = time.time()
    global path_find_time_consume
    path_find_time_consume += ed - st

    paths_list.sort(key=sort_key)

    for path,weight,use_contig in paths_list:
        # below: create copy of path with each node as rc version
        # use as unique representation of a path and rc of its whole
        unoriented_sorted_path_str = get_unoriented_sorted_str(path)

        # here we avoid considering cyclic rotations of identical paths
        # by sorting their string representations
        # and comparing against the set already stored
        if unoriented_sorted_path_str not in unq_sorted_paths:
            unq_sorted_paths.add(unoriented_sorted_path_str)
            paths.append(tuple((path,weight,use_contig)))

    return paths


def get_high_mass_shortest_path(node,G,path_dict,SEQS,use_scores,use_genes):
    """ Return the shortest circular path back to node
    """
    # rc_node_hit used to determine whether use node or rc_node(node) to find plasmid gene hit cycle

    # TODO: potentially add check for unique paths so that don't check same cycle
    # twice if there are two potential plasmid nodes in it

    for e in G.edges():
        if use_genes and G.nodes[e[1]]['gene'] == True:
            G.add_edge(e[0], e[1], cost = 0.0)
        elif use_scores == True:
            G.add_edge(e[0], e[1], cost = (1.-(G.nodes[e[1]]['score']))/get_spades_base_mass(G, e[1]))
        else:
            G.add_edge(e[0], e[1], cost = (1./get_spades_base_mass(G, e[1])))

    shortest_score = float("inf")
    path = None
    shortest_path = None
    use_contig = False
    
    st = time.time()
    for pred in G.predecessors(node):
        try:
            length, path, through_contig = dijkstra_path(G,path_dict,SEQS,source=node,target=pred, weight='cost')
            if length < shortest_score:
                shortest_path = tuple(path)
                shortest_score = length
                use_contig = through_contig
            # logger.info(f"special path:  {str(folded_path)}")
        except nx.exception.NetworkXNoPath:
            continue
    ed = time.time()
    global path_find_time_consume
    path_find_time_consume += ed - st
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


def count_selfloop_mates(node,bamfile,G):
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

####### Updated - exclude path with node that mostly has mates outside of path
def is_good_cyc(path, G, bamfile):
    """ check all non-repeat nodes only have mates
        mapping to contigs in the cycle
    """
    # TODO we may need use this construction in final step
    # return True

    sing_nodes = set()
    for node in path:
      if node[-1] == "'": node = node[:-1]
      sing_nodes.add(node)
    if len(sing_nodes)==0: return True

    non_path_dominated_nodes = 0

    for nd in sing_nodes:
        mate_tigs = get_contigs_of_mates(nd, bamfile, G)
        # NOTE: ^ this only gets mates that are reachable from nd in G

        # logger.info("\tNode: %s" % nd)
        # logger.info("\t\tMates: %s" % ", ".join(mate_tigs))

        # need to check against F and R versions of path nodes
        path_rc = [rc_node(x) for x in path]
        num_mates_in_path = sum([1 for x in mate_tigs if (x in path or x in path_rc)])
        num_mates_not_in_path = len(mate_tigs)-num_mates_in_path
        if num_mates_in_path < num_mates_not_in_path:
       ########### if len(mate_tigs)>1 and num_mates_in_path < num_mates_not_in_path:
            non_path_dominated_nodes += 1
    # if float(non_path_dominated_nodes)/float(len(sing_nodes)) >= PARAMS.GOOD_CYC_DOMINATED_THRESH:
    if non_path_dominated_nodes == len(sing_nodes):


        logger.info(f"Too many nodes with majority of mates not on path: {non_path_dominated_nodes}, {len(sing_nodes)}")
        return False
    else: return True
    
#########################
def process_component(COMP, G, original_comp, max_k, min_length, max_CV, SEQS, bamfile, pool, path_dict,node_to_contig,contigs_path_name_dict, use_scores=False, use_genes=False, num_procs=1):
    """ run recycler for a single component of the graph
        use multiprocessing to process components in parallel
    """
        ###############MOVED FROM OUTER CODE ON WHOLE G
        # may be redundant
    if use_scores: remove_hi_confidence_chromosome(COMP) ##################################

    # initialize shortest path set considered
    path_count = 0
    seen_unoriented_paths = set([])
    paths_set = set([]) #the set of paths found

    # looking for paths starting from the nodes annotated with plasmid genes
    if use_genes:
        plasmid_gene_nodes = get_plasmid_gene_nodes(COMP)
        potential_plasmid_mass_tuples = [(get_spades_base_mass(COMP,nd),nd) for nd in plasmid_gene_nodes]
        potential_plasmid_mass_tuples.sort(key = lambda n: n[0]) # sorted by coverage * length
        while potential_plasmid_mass_tuples: # could be removing other nodes from the list
            top_node = potential_plasmid_mass_tuples.pop() # highest mass node
            top_node_name = top_node[1]
            logger.info(f"plamsid gene hit node: {top_node_name}")
            # 初始化rc_node的相关信息
            rc_length, rc_path,rc_use_contig,rc_path_CV = float("inf"), None, False,float("inf")
            length, path,use_contig = get_high_mass_shortest_path(top_node_name,COMP,path_dict,SEQS,use_scores,use_genes) #######
            path_CV = float("inf")
            if path is not None:
                path_CV = get_wgtd_path_coverage_CV(path,G,SEQS,max_k_val=max_k)
                logger.info(f"plasmid gene path: {path}")
                logger.info(f"CV: {path_CV}, Good: {is_good_cyc(path,G,bamfile)}, weight: {length}, use_contig: {use_contig}")
            if rc_node(top_node_name) in plasmid_gene_nodes:
                rc_length, rc_path,rc_use_contig = get_high_mass_shortest_path(rc_node(top_node_name),COMP,path_dict,SEQS,use_scores,use_genes)
                if rc_path is not None:
                    rc_path_CV = get_wgtd_path_coverage_CV(rc_path,G,SEQS,max_k_val=max_k)
                    logger.info(f"RC plasmid gene path: {rc_path}") 
                    logger.info(f"CV: {rc_path_CV}, Good: {is_good_cyc(rc_path,G,bamfile)}, weight: {rc_length}, use_contig:{rc_use_contig}")
            
            if path is None and rc_path is None: continue
            # 由于正向序列或反向互补序列第一次被添加后就会删除另一个，因此需要通过逻辑判断确定最终的结果：
            # 如果反向互补序列的最短路径更短，并且符合添加条件，则使用反向互补序列，否则使用正向序列
            if path is None:
                path = rc_path
            elif rc_path is not None:
                if meet_criterion(path, G, SEQS,max_k, max_CV,bamfile):
                    if meet_criterion(rc_path, G, SEQS,max_k, max_CV,bamfile):
                        if rc_use_contig > use_contig or (rc_use_contig == use_contig and rc_path_CV < path_CV):
                            path = rc_path
                elif meet_criterion(rc_path, G, SEQS,max_k, max_CV,bamfile):
                    path = rc_path         
                    
            if meet_criterion(path, G, SEQS,max_k, max_CV,bamfile):
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
                paths_set.add((path,before_cov))
            else:
                logger.info("Did not add plasmid gene path: %s" % (str(path)))

                
                
            
        # then look for circular paths that start from hi confidence plasmid nodes
    if use_scores:
        potential_plasmid_nodes = get_hi_conf_plasmids(COMP)
        potential_plasmid_mass_tuples = [(get_spades_base_mass(COMP,nd),nd) for nd in potential_plasmid_nodes]
        potential_plasmid_mass_tuples.sort(key = lambda n: n[0])
        while potential_plasmid_mass_tuples: # could be removing other nodes from the list
            top_node = potential_plasmid_mass_tuples.pop() # highest mass node
            top_node_name = top_node[1]
            logger.info(f"Hi conf node: {top_node_name}")
            rc_length, rc_path,rc_use_contig,rc_path_CV = float("inf"), None, False,float("inf")
            length, path,use_contig = get_high_mass_shortest_path(top_node_name,COMP,path_dict,SEQS,use_scores,use_genes) #######
            path_CV = float("inf")
            if path is not None:
                path_CV = get_wgtd_path_coverage_CV(path,G,SEQS,max_k_val=max_k)
                logger.info(f"Hi conf path: {path}") 
                logger.info(f"CV: {path_CV}, Good: {is_good_cyc(path,G,bamfile)}, weight: {length}, use_contig: {use_contig}")
                
            if rc_node(top_node_name) in potential_plasmid_nodes:
                rc_length, rc_path,rc_use_contig = get_high_mass_shortest_path(rc_node(top_node_name),COMP,path_dict,SEQS,use_scores,use_genes)
                if rc_path is not None:
                    rc_path_CV = get_wgtd_path_coverage_CV(rc_path,G,SEQS,max_k_val=max_k)
                    logger.info(f"RC Hi conf path: {rc_path}")
                    logger.info(f"CV: {rc_path_CV}, Good: {is_good_cyc(rc_path,G,bamfile)}, weight: {rc_length}, use_contig: {rc_use_contig}")
            
            if path is None and rc_path is None: continue
            # 由于正向序列或反向互补序列第一次被添加后就会删除另一个，因此需要通过逻辑判断确定最终的结果：
            # 如果反向互补序列的最短路径更短，并且符合添加条件，则使用反向互补序列，否则使用正向序列
            if path is None:
                path = rc_path
            elif rc_path is not None and path is not None:
                if meet_criterion(path, G, SEQS,max_k, max_CV,bamfile):
                    if meet_criterion(rc_path, G, SEQS,max_k, max_CV,bamfile):
                        if rc_use_contig > use_contig or (rc_use_contig == use_contig and rc_path_CV < path_CV):
                            path = rc_path
                elif meet_criterion(rc_path, G, SEQS,max_k, max_CV,bamfile):
                    path = rc_path         
                    
            if meet_criterion(path, G, SEQS,max_k, max_CV,bamfile):
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
                paths_set.add((path,before_cov))
            else:
                logger.info("Did not add Hi conf path: %s" % (str(path)))


        # 3rd step. Run Recycler algorithm that looks for circular high mass shortest
        # paths and accept them as plasmid predictions if the coverages and mate pairs
        # match the required thresholds
#######################################################################################
#######################################################################################


    paths = enum_high_mass_shortest_paths(COMP, pool,path_dict ,SEQS,use_scores,use_genes,seen_unoriented_paths)
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
        for p,length,use_contig in paths:
            # if len(get_seq_from_path(p, SEQS, max_k_val=max_k)) < min_length:
            if get_total_len_from_path(p,max_k_val=max_k,cycle=True) < min_length:
                seen_unoriented_paths.add(get_unoriented_sorted_str(p))
                # logger.info("Num seen paths: %d" % (len(seen_unoriented_paths)))
                continue
            path_tuples.append((get_wgtd_path_coverage_CV(p,G,SEQS,max_k_val=max_k), p,length,use_contig))
            

        logger.info("Num path tuples: %d" % (len(path_tuples)))
        if(len(path_tuples)==0): break

        # 使用了contig path的优先，低CV作为第二排序依据
        path_tuples.sort(key=lambda item: (not item[3], item[0]))

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
                    is_good_cyc(curr_path,G,bamfile)):

                    logger.info("Added path %s" % ", ".join(curr_path))
                    logger.info("\tCV: %4f" % curr_path_CV)
                    seen_unoriented_paths.add(get_unoriented_sorted_str(curr_path))
                    #before_cov, _ = get_path_mean_std(curr_path, COMP, SEQS, max_k)
                    before_cov, _ = get_path_mean_std(curr_path, G, SEQS, max_k)
                    covs = update_path_coverage_vals(curr_path, G, SEQS, max_k)
                    update_path_with_covs(curr_path, COMP, covs)
                    path_count += 1
                    paths_set.add((tuple(curr_path), before_cov))
                    break

                else:
                    logger.info("Did not add path: %s" % (", ".join(curr_path)))
                    logger.info("\tCV: %4f" % curr_path_CV)
                    if curr_path_CV > max_CV:
                        break # sorted by CV
                    else: # not good mate pairs
                        seen_unoriented_paths.add(get_unoriented_sorted_str(curr_path))

        # recalculate paths on the component
        print(str(len(COMP.nodes())) + " nodes remain in component")
        logger.info("Remaining nodes: %d" % (len(COMP.nodes())))
        paths = enum_high_mass_shortest_paths(COMP, pool,path_dict,SEQS, use_scores,use_genes,seen_unoriented_paths)

    # #end while
    st = time.time()
    path_set_before_merge = paths_set.copy()
    logger.info("start merge...")
    merged_paths_set = merge_cycle(paths_set,SEQS,max_k,original_comp,bamfile,node_to_contig,contigs_path_name_dict)
    logger.info("merge finished:")
    # merged_paths_set = paths_set
    ed = time.time()
    merge_time = ed - st
    return merged_paths_set, merge_time, path_set_before_merge

            
def merge_cycle(paths_set:set,SEQS,max_k,G,bamfile,node_to_contig,contigs_path_name_dict):
    #记录原有set规模，如果有合并，则规模改变
    prev_len = -1
    #由于此component为强连通分支，如果正向和反向互补序列均在图中，此时选择的cycle已经是两者中更好的，如果反向互补序列不在此连通分支中，则它们理应不连通，也无需考虑
    score_dict = get_score_from_set(paths_set,SEQS)
    while(prev_len != len(paths_set)):
        prev_len = len(paths_set)
        # path_set中的元素为 (path, cov)
        for cur_path,cur_cov in sorted(paths_set,key=lambda x:x[1],reverse=True):
           
            cur_set = set(cur_path)
            cur_score = score_dict[(cur_path,cur_cov)]

            for other_path, other_cov in sorted(paths_set,key=lambda x:x[1],reverse=True):
                if other_path == cur_path:
                    continue
                original_other_path = other_path
                other_set = set(other_path)
                # rc_other_set = set(rc_path(other_path))
                intersection = cur_set.intersection(other_set)
                # if rc_other_set.issubset(SEQS):
                #     # TODO we need to take RC into consideration.
                #     rc_intersction = cur_set.intersection(rc_other_set)
                #     if len(rc_intersction) > len(intersection):
                #         other_path = rc_path(other_path)
                #         intersection = rc_intersction

                if len(intersection)>0:
                    logger.info(f"cur_pth: {cur_path}\nother_path: {other_path}\nitersection: {intersection}")
                    to_merge, to_double_merge =  is_similar(cur_path,other_path,cur_cov,other_cov,SEQS,max_k)
                    merged_path = None
                    # 依据 contig path info, contig 只考虑单倍cov
                    if to_merge:
                        for nd in sorted(intersection, key=lambda nd: get_cov_from_spades_name(nd), reverse=True):
                            if nd in node_to_contig.keys():
                                merged_path = merge_path_through_contig_path(cur_path,other_path,nd, node_to_contig,contigs_path_name_dict)
                                if merged_path is not None:
                                    break

                    if merged_path is not None:
                        merged_score = get_score_from_path(merged_path,SEQS)
                        merged_cv = get_wgtd_path_coverage_CV(merged_path,G,SEQS,max_k_val=max_k)
                        logger.info(f"curpath: {cur_path}, cov: {cur_cov}, score: {score_dict[(cur_path,cur_cov)]}")
                        logger.info(f"\totherpath: {other_path}, cov: {other_cov}, score: {score_dict[(other_path,other_cov)]}")
                        merged_cov , merged_std = get_path_mean_std(merged_path, G, SEQS, max_k,discount=True)
                        logger.info(f"merged_path_through_contig: {merged_path}, cov: {merged_cov}, CV: {merged_cv}, score: {merged_score}")
                        if merged_score >= 0.5 and merged_cv <= 0.5:
                            paths_set.add((merged_path,merged_cov))
                            score_dict[(merged_path,merged_cov)] = merged_score
                            # score_dict[(rc_path(merged_path),merged_cov)] = get_score_from_path(rc_path(merged_path),SEQS)
                            #删除原有path
                            paths_set.remove((cur_path,cur_cov))
                            paths_set.remove((original_other_path,other_cov))
                            score_dict.pop((cur_path,cur_cov))
                            score_dict.pop((original_other_path,other_cov))
                            # score_dict.pop((rc_path(cur_path),cur_cov))
                            # score_dict.pop((rc_path(original_other_path),other_cov))
                           
                            break
                        else:
                            logger.info(f"bad contig merge: merged_score: {merged_score}, merged_cv:{merged_cv}, skip")

                    
                    if to_merge or to_double_merge:
                        # Supposing that a cycle is a confident contig,namely,we randamly select a enterpoint and insert the whole path
                        enterpoint = sorted(intersection, key=lambda nd: get_cov_from_spades_name(nd), reverse=True)[0]
                        logger.info(f" enterpoint is {enterpoint}")
                        logger.info(f" attempt to merge :\n{cur_path}\n{other_path}")
                        other_score = score_dict[(other_path,other_cov)]
                        merged_path = merge_paths(cur_path,other_path, enterpoint, to_double_merge)
                        merged_score = get_score_from_path(merged_path,SEQS)
                        merged_cv = get_wgtd_path_coverage_CV(merged_path,G,SEQS,max_k_val=max_k)
                        # TODO 判断融合前后的score是否提升了，没提升则不合并
                        # and get_wgtd_path_coverage_CV(p,G,SEQS,max_k_val=max_k) < 0.5
                        # potential bugs: c1_score : 0.9, c2_score :0.001 --> merged_score 0.5
                        # c1 合并c2 不成功，c2 合并c1 会成功
                        if merged_score >= cur_score  and merged_score >= other_score and merged_cv <= 0.5:
                        # if merged_score >= cur_score*0.95 and merged_score >= other_score*0.95  and merged_cv < 0.8:

                            # 合并
                            logger.info(f"curpath: {cur_path} \ncov: {cur_cov}, score: {cur_score}")
                            logger.info(f"\totherpath: {other_path} \ncov: {other_cov}, score: {other_score}")
                            merged_cov , merged_std = get_path_mean_std(merged_path, G, SEQS, max_k,discount=True)
                            logger.info(f"merged_path: {merged_path} \ncov: {merged_cov}, CV: {merged_cv}, score: {merged_score}")
                            paths_set.add((merged_path,merged_cov))
                            score_dict[(merged_path,merged_cov)] = merged_score
                            # score_dict[(rc_path(merged_path),merged_cov)] = get_score_from_path(rc_path(merged_path),SEQS)
                            #删除原有path
                            paths_set.remove((cur_path,cur_cov))
                            paths_set.remove((original_other_path,other_cov))
                            score_dict.pop((cur_path,cur_cov))
                            score_dict.pop((original_other_path,other_cov))
                            # score_dict.pop((rc_path(cur_path),cur_cov))
                            # score_dict.pop((rc_path(original_other_path),other_cov))
                            #直接重新进行，防止在循环中修改set可能导致的bug
                            break

                        else:
                            logger.info(f"bad merge: cur_score: {cur_score}, other_score: {other_score} --> merged_score: {merged_score}, merged_cv:{merged_cv}, skip")
                    else:
                        logger.info(f"not similar, skip")
            
            if prev_len != len(paths_set):
                #直接重新进行，防止在循环中修改set可能导致的bug
                break
            
    good_paths_set = set()

    for path, cov in paths_set:
        logger.info(f"Path: {path}")
        if is_good_cyc(path,G,bamfile):
            good_paths_set.add((path,cov))   
    return  good_paths_set

def get_score_from_set(paths_set: set, SEQs):
    score_dict = {}
    for path,cov in paths_set:
        score = get_score_from_path(path,SEQs)
        score_dict[(path,cov)] = score
        # rc_p = rc_path(path)
        # if all(nd in SEQs for nd in rc_p):
        #     rc_score = get_score_from_path(rc_p,SEQs)
        #     score_dict[(rc_p,cov)] = rc_score
    return score_dict

def merge_paths(path1, path2,enterpoint,is_double):
    
    # 将短序列插入到长序列中
    if len(path1) < len(path2):
        path1 , path2 = path2, path1

    # 找到路径中交集元素的位置
    index1 = path1.index(enterpoint)
    index2 = path2.index(enterpoint)

    # 合并路径：合并顺序按照 path1 + path2 的顺序，并去除重复节点
    if is_double:
        logger.info(f"double merge")
        merged_path = path1[:index1] + path2[index2:] + path2[:index2] + path2[index2:] + path2[:index2] + path1[index1:]
    else:
        merged_path = path1[:index1] + path2[index2:] + path2[:index2] + path1[index1:]
    return merged_path
    
def is_similar(cur_path:list,other_path:list,cur_set_mean_cov,other_set_mean_cov,SEQS,max_k, cutoff=0.15, factor=0.8):
    # # 重新计算不包含重复节点的coverage
    # cur_set = set(cur_path)
    # other_set= set(other_path)
    # intersection = cur_set.intersection(other_set)
    # cur_set-=intersection
    # other_set-=intersection
    # cur_set_mean_cov = get_mean_cov_from_set(cur_set)
    # other_set_mean_cov = get_mean_cov_from_set(other_set)


    # 考虑小环可能经过多次，所以小环的cov会是大环的倍数，合格后的cov应该选择大环的mean_cov
    cur_set_mean_cov = float(cur_set_mean_cov)
    other_set_mean_cov = float(other_set_mean_cov)
    cur_path_len = len(get_seq_from_path(cur_path, SEQS, max_k_val=max_k))
    other_path_len = len(get_seq_from_path(other_path, SEQS, max_k_val=max_k))
    shorter_path_cov, longger_path_cov = 0, 0 
    if cur_path_len < other_path_len:
        shorter_path_cov, longger_path_cov = cur_set_mean_cov,other_set_mean_cov
        longger_path_len, shorter_path_len = other_path_len, cur_path_len
    else:
        shorter_path_cov,longger_path_cov  = other_set_mean_cov, cur_set_mean_cov
        shorter_path_len, longger_path_len = other_path_len, cur_path_len
    # 两环mean coverage相近
    similar_difference = get_relative_difference(shorter_path_cov,longger_path_cov)

    # 小环的mean coverage是大环的倍数(这里取2)
    double_similar_difference = get_relative_difference(shorter_path_cov , longger_path_cov * 2)
    logger.info(f"cur_path_len: {cur_path_len}, other_path_len: {other_path_len}")
    logger.info(f"similar_diff: {similar_difference}, double_similar_diff: {double_similar_difference}")
    # TODO 当cov是倍数关系时，最终环的cov选择长环的，否则选两者的均值
    return  (similar_difference < cutoff)  , (double_similar_difference < cutoff and shorter_path_len/longger_path_len < factor)

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

def get_score_from_path(path,SEQS,num_procs=1):
    seqs = get_seq_from_path(path,SEQS)
    c = plasclass.plasclass(num_procs)
    prob = c.classify(seqs)
    return prob

def get_contig_path(path_file,id_dict,SEQS,G,contig_path_file,score_out_file,min_contig_path_len=4,max_k=77,num_procs=1):
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

                    if len(path_list) < min_contig_path_len:
                        continue

                    # 自定义的contig path的唯一标识
                    unq_contig_name = contig_name+"_"+str(count)
                        

                    # add reversed path for bidireactional dijistra
                    contigs_path_name_dict[unq_contig_name] = path_list
                    reversed_path_list = list(reversed(path_list))
                    contigs_path_name_dict['R'+unq_contig_name] = reversed_path_list
                    count+=1

                    # store contig path 
                    # the judgement may be redundant
                    if all(node in node_set for node in path_list ):
                        contigs_path_dict[0].setdefault(path_list[0],[])
                        contigs_path_dict[1].setdefault(reversed_path_list[0],[])

                        # dont store the first node which stored in key
                        # key:  first node of contig path
                        # value: (contig path, contig path name, score of contig path, conitg path fragment end with key(here is Null))
                        contigs_path_dict[0][path_list[0]].append((path_list[1:],unq_contig_name,0,None))
                        contigs_path_dict[1][reversed_path_list[0]].append((reversed_path_list[1:],'R'+unq_contig_name,0,None))
                            
        logger.info(f"preprocess contig path score...")
        print("Preprocessing contig path score")
        start_time = time.time()

        # 将所有的contig path写入文件并批量预测得分
        for k in (0,1):
            for node, info in contigs_path_dict[k].items():
                for path, name,_,_ in info:
                    if len(path) >= min_contig_path_len-1:
                        o.write(f">{name}\n")
                        new_path = [node]+path[:]
                        o.write(f"{get_seq_from_path(new_path,SEQS,max_k,False)}\n")  
                    else:
                        raise ValueError(f"path of {node} have illeagel length")
    classify(contig_path_file, score_out_file, num_procs)

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
            for idx, (path, name, s, pre) in enumerate(info_list):
                if len(path) >= min_contig_path_len - 1:
                    score = scores_dict.get(name, None)
                    if score is not None and score >= 0.5:
                        info_list[idx] = (path, name, score, None)  # 更新整个元组
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


def get_contig_path_cost(path,G,seqs,score,discount=True,max_k_val=77):
    for n in path:
        if G.nodes[n]['gene'] == True:
            return 0.0
    try:
        covs = np.array(get_path_covs(path,G,discount))
        wgts = np.array([(get_length_from_spades_name(n)-max_k_val) for n in path])
        # TODO can be obtained by name and number of nodes
        # seq  = get_seq_from_path(path, seqs, max_k_val, cycle=False)
        tot_len = get_total_len_from_path(path,max_k_val=max_k_val)

        wgts = np.multiply(wgts, 1./tot_len)
        mean = np.average(covs, weights = wgts)
    except ZeroDivisionError as e:
        print(f"path:{path}")
        sys.exit(f"contig path length is zero")

    return (1-score)/mean * tot_len
    
    
def dijkstra_path(G, path_dict:dict,SEQS, source, target, weight='weight',bidirectional=False):


    (length, path,use_contig) = single_source_dijkstra(G, path_dict,SEQS ,source, target=target,weight=weight,bidirectional=bidirectional)
    return length,path,use_contig

def single_source_dijkstra(G, path_dict,SEQS,source, target=None, cutoff=None,
                           weight='weight',bidirectional=False):
    

    return multi_source_dijkstra(G, path_dict,SEQS,[source], cutoff=cutoff, target=target,
                                 weight=weight,bidirectional=bidirectional)

def multi_source_dijkstra(G, path_dict,SEQS,sources, target=None, cutoff=None,
                          weight='weight',bidirectional=False):

    use_contig = False
    if not sources:
        raise ValueError('sources must not be empty')
    if target in sources:
        return (0, [target],use_contig)
    weight_lambda = _weight_function(G, weight)
    paths = {source: [source] for source in sources}  # dictionary of paths
    # folded_paths = {source: [source] for source in sources} 
    G_succ = G._succ if G.is_directed() else G._adj
    node_set =set(G_succ.keys())
    min_len = float("inf")
    final_path = None
    if target is None:
        return (min_len, paths,use_contig)
   
    # 单向dijstra只利用path_dict中存储的正向的contig path，即path_dict[0]
    # sources[0] 就是source, 因为传入的是一个数组
    if sources[0] in path_dict[0].keys():
        # logger.info(f"source: {sources[0]}")
        for record in path_dict[0][sources[0]]:
            try:
               
                path,contig_name,score,pre_contig = record
                cur_path = {path[-1] : [path[-1]]}
                # logger.info(f"path dict has {contig_name}")
                finaldist, finalpath = float("inf"),None
                # 当前节点为contig path的首节点，直接找contig path的尾节点到当前节点的前驱节点的最短路径即可
                if pre_contig is None:
                    if all(node in node_set for node in path):
                        
                        if bidirectional is False:
                            finaldist, finalpath,_ = _dijkstra_multisource(G, path_dict,SEQS,[path[-1]], weight_lambda, paths=cur_path,
                                        cutoff=cutoff, target=target)
                        else:
                            finaldist, finalpath,_ = bidirectional_dijkstra(G, path_dict,SEQS,path[-1], target,weight)

                        finaldist+= get_contig_path_cost([sources[0]]+path,G,SEQS,score,discount=True,max_k_val=77)
                        # add first contig to final path
                        finalpath = [sources[0]]+ path + finalpath[1:]
                        # logger.info(f"use contig {contig_name}")

                    else:
                        # logger.info(f"{contig_name} not in G: {', '.join(str(node) for node in path if node not in node_set)}")
                        continue
                # source is in contig path
                # 当前节点在contig path中，直接找contig path的尾节点到cotig path的首节点的最短路径即可
                else:
                    if all(node in node_set for node in path) and all(node in node_set for node in pre_contig):
                        if bidirectional is False:
                            finaldist, finalpath, _ = _dijkstra_multisource(G, path_dict,SEQS,[path[-1]], weight_lambda, paths=cur_path,
                                        cutoff=cutoff, target=pre_contig[0])
                        else:
                            finaldist, finalpath,_ = bidirectional_dijkstra(G, path_dict,SEQS,path[-1], pre_contig[0],weight)

                        finaldist = get_contig_path_cost(pre_contig+[sources[0]]+path,G,SEQS,score,discount=True,max_k_val=77) + finaldist
                        finalpath = [sources[0]]+path + finalpath[1:] + pre_contig[1:]
                        # logger.info(f"use contig {contig_name}")

                    else:
                        # logger.info(f"{contig_name} not in G: {', '.join(str(node) for node in pre_contig+path if node not in node_set)}")
                        continue
                # logger.info(f"contig path: {sources[0]}--->{target}: {finalpath}, weight: {finaldist}")
                if finaldist < min_len:
                    min_len = finaldist
                    final_path = finalpath
                    # logger.info(f"final add contig path: {sources[0]}--->{target}: {finalpath}, weight: {finaldist}")
                    use_contig = True
            except KeyError:
                # logger.info(f"has not path from {sources[0]} to {target} through {contig_name}")
                continue
    # can't find path through contig
    # 如果通过contig找不到最短环，则直接通过双向dijistra找
    if final_path is None:
        finaldist, finalpath, through_contig= bidirectional_dijkstra(G, path_dict,SEQS,sources[0], target,weight)
        # finaldist, finalpath, through_contig= _dijkstra_multisource(G, path_dict,SEQS,sources, weight_lambda, paths=paths,
        #                          cutoff=cutoff, target=target)
        try:
            # logger.info(f"simple path: {sources[0]}--->{target}: {finalpath}, weight: {finaldist}")
            min_len = finaldist
            final_path = finalpath
            use_contig = through_contig
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


def _dijkstra_multisource(G, path_dict:dict,SEQS,sources, weight,pred=None, paths=None,
                          cutoff=None, target=None):
    # TODO 如果被一个以repeat起始的contig path 引导进入了一个 subloop, 出来时再度遇到repeat会导致算法直接结束
    G_succ = G._succ if G.is_directed() else G._adj
    node_set = set(G_succ.keys())
    push = heappush
    pop = heappop
    dist = {}  # dictionary of final distances
    seen = {}
    # fringe is heapq with 3-tuples (distance,c,node)
    # use the count c to avoid comparing nodes (may not be able to)
    c = count()
    use_contig = {}
    fringe = []
    for source in sources:
        if source not in G:
            raise nx.NodeNotFound("Source {} not in G".format(source))
        seen[source] = 0
        use_contig[source] = False
        push(fringe, (0, next(c), source))
    while fringe:
        through_contig = False
        (d, _, v) = pop(fringe)
        if v in dist:
            continue  # already searched this node.
        dist[v] = d
        if v == target:
            break
        if v in path_dict[0].keys():
            for record in path_dict[0][v]:
                path,contig_name,score,pre_contig = record
                if pre_contig is not None:
                    continue
                if all(node in node_set for node in path):
                    cost = get_contig_path_cost(path,G,SEQS,score)
                    if cost < 0 :
                        raise ValueError(f"contig path has negative weight {path,cost}")
                    
                    vu_dist = dist[v]+cost
                    if cutoff is not None:
                        if vu_dist > cutoff:
                            continue
                    u = path[-1]
                    if u in dist:
                        if vu_dist < dist[u]:
                            raise ValueError('Contradictory paths found:','negative weights?')
                    elif u not in seen or vu_dist < seen[u]:
                        through_contig = True
                        seen[u] = vu_dist
                        push(fringe, (vu_dist, next(c), u))
                        if paths is not None:
                            paths[u] = paths[v] + path[:]
                            use_contig[u]  = use_contig[v] | through_contig
                            # folded_paths[u] = folded_paths[v] + [record[1]]
                            
                        if pred is not None:
                            pred[u] = [v]
                    elif vu_dist == seen[u]:
                        if pred is not None:
                            pred[u].append(v)
        if through_contig == False:
            for u, e in G_succ[v].items():
                cost = weight(v, u, e)
                if cost is None:
                    continue
                vu_dist = dist[v] + cost
                if cutoff is not None:
                    if vu_dist > cutoff:
                        continue
                if u in dist:
                    if vu_dist < dist[u]:
                        raise ValueError('Contradictory paths found:',
    'negative weights?')
                elif u not in seen or vu_dist < seen[u]:
                    seen[u] = vu_dist
                    push(fringe, (vu_dist, next(c), u))
                    if paths is not None:
                        paths[u] = paths[v] + [u]
                        use_contig[u]  = use_contig[v] | through_contig
                        # folded_paths[u] = folded_paths[v] + [u]
                    if pred is not None:
                        pred[u] = [v]
                elif vu_dist == seen[u]:
                    if pred is not None:
                        pred[u].append(v)
    try:
        finaldist=dist[target]
        finalpath=paths[target]
        Use_contig=use_contig[target]
    except KeyError:
        raise nx.NetworkXNoPath("No path between %s and %s." % (source[0], target))

    return finaldist,finalpath,Use_contig

def bidirectional_dijkstra(G,  path_dict,SEQS, source, target, weight='weight'):

    if source not in G or target not in G:
        msg = 'Either source {} or target {} is not in G'
        raise nx.NodeNotFound(msg.format(source, target))

    if source == target:
        return (0, [source],False)
    push = heappush
    pop = heappop
    # Init:  [Forward, Backward]
    dists = [{}, {}]   # dictionary of final distances
    paths = [{source: [source]}, {target: [target]}]  # dictionary of paths
    use_contig = [{source: False},{target: False}]
    fringe = [[], []]  # heap of (distance, node) for choosing node to expand
    seen = [{source: 0}, {target: 0}]  # dict of distances to seen nodes
    c = count()
    # initialize fringe heap
    push(fringe[0], (0, next(c), source))
    push(fringe[1], (0, next(c), target))
    # neighs for extracting correct neighbor information
    if G.is_directed():
        neighs = [G.successors, G.predecessors]
    else:
        neighs = [G.neighbors, G.neighbors]
    # variables to hold shortest discovered path
    finaldist = float("inf")
    finalpath = []
    dir = 1
    node_set = set(G.nodes())
    while fringe[0] and fringe[1]:
        through_contig = False
        # choose direction
        # dir == 0 is forward direction and dir == 1 is back
        dir = 1 - dir
        # extract closest to expand
        (dist, _, v) = pop(fringe[dir])
        if v in dists[dir]:
            # Shortest path to v has already been found
            continue
        # update distance
        dists[dir][v] = dist  # equal to seen[dir][v]
        if v in dists[1 - dir]:
            # if we have scanned v in both directions we are done
            # we have now discovered the shortest path
            return (finaldist, finalpath,use_contig[dir][v] | use_contig[1-dir][v])
        if v in path_dict[dir].keys():
            for record in path_dict[dir][v]:
                path, contig_name,score,pre_contig = record
                if pre_contig is not None:
                    continue
                if all(node in node_set for node in path):
                    # TODO max_k 需要传入
                    
                    costs = get_contig_path_cost(path,G,SEQS,score)
                    w = path[-1]
                    
                    vwLength = dists[dir][v] + costs
                    if w in dists[dir]:
                        if vwLength < dists[dir][w]:
                            raise ValueError(
                                "Contradictory paths found: negative weights?")
                    elif w not in seen[dir] or vwLength < seen[dir][w]:
                        through_contig = True
                        seen[dir][w] = vwLength
                        push(fringe[dir], (vwLength, next(c), w))
                        paths[dir][w] = paths[dir][v] + path[:]
                        use_contig[dir][w] = use_contig[dir][v] | through_contig
                        if w in seen[0] and w in seen[1]:
                            # see if this path is better than than the already
                            # discovered shortest path
                            totaldist = seen[0][w] + seen[1][w]
                            if finalpath == [] or finaldist > totaldist:
                                finaldist = totaldist
                                revpath = paths[1][w][:]
                                revpath.reverse()
                                finalpath = paths[0][w] + revpath[1:]

        if through_contig == False:
            for w in neighs[dir](v):
                if(dir == 0):  # forward
                    if G.is_multigraph():
                        minweight = min((dd.get(weight, 1)
                                        for k, dd in G[v][w].items()))
                    else:
                        minweight = G[v][w].get(weight, 1)
                    vwLength = dists[dir][v] + minweight  # G[v][w].get(weight,1)
                else:  # back, must remember to change v,w->w,v
                    if G.is_multigraph():
                        minweight = min((dd.get(weight, 1)
                                        for k, dd in G[w][v].items()))
                    else:
                        minweight = G[w][v].get(weight, 1)
                    vwLength = dists[dir][v] + minweight  # G[w][v].get(weight,1)

                if w in dists[dir]:
                    if vwLength < dists[dir][w]:
                        raise ValueError(
                            "Contradictory paths found: negative weights?")
                elif w not in seen[dir] or vwLength < seen[dir][w]:
                    # relaxing
                    seen[dir][w] = vwLength
                    push(fringe[dir], (vwLength, next(c), w))
                    paths[dir][w] = paths[dir][v] + [w]
                    use_contig[dir][w] = use_contig[dir][v] | through_contig
                    if w in seen[0] and w in seen[1]:
                        # see if this path is better than than the already
                        # discovered shortest path
                        totaldist = seen[0][w] + seen[1][w]
                        if finalpath == [] or finaldist > totaldist:
                            finaldist = totaldist
                            revpath = paths[1][w][:]
                            revpath.reverse()
                            finalpath = paths[0][w] + revpath[1:]
    raise nx.NetworkXNoPath("No path between %s and %s." % (source, target))

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


def add_contig_to_path_dict(G,scores_dict, path_dict, contigs_path_name_dict,node_to_contig, use_genes=True, use_scores=True):

    # 对于处在contig path中间的 gene hit node 或 hi conf node, 单独新增一些记录：如 3 为 hi conf node
    # 之前存储过cotig path : {key: 1, value: ([2 3 4 5 6],contig_name, score, None) }   
    # 现在新增记录: {key: 3, value:([4 ,5 ,6], contig_name, score, [1, 2])} 
    # 新纪录用于将通过hi conf node的 contig path 完整利用起来，如以3为起点寻找最短路径，则可以直接找 6 - > 1的最短路径，再拼上整条路径
    plasmid_nodes = set()

    # Determine plasmid nodes based on genes and scores
    if use_genes:
        gene_hit_nodes = get_plasmid_gene_nodes(G)  
        plasmid_nodes.update(gene_hit_nodes)
    if use_scores:
        hi_conf_nodes = get_hi_conf_plasmids(G)
        plasmid_nodes.update(hi_conf_nodes)

    for node in plasmid_nodes:
        if node not in node_to_contig.keys():
            logger.info(f"hi conf node: {node} not in contig path")
            #node not in contig path
            continue
        for contig_id in node_to_contig[node]:
            # 对于gene hit node 和 hi conf node， 使用单项dijstra,因此无需考虑反向序列
            if contig_id.startswith('R'):
                continue
            # 通过contig name 获取 该node所在的path
            path = contigs_path_name_dict[contig_id]
            if node in path:  # Ensure the node is in the path before indexing
                start_index = path.index(node)

                # 如果该节点已经位于开头或者末尾，则无需新增path
                # TODO 尾部可能需要单独处理
                if start_index == len(path)-1 or start_index == 0 :
                    continue
                contig_path1 = path[start_index:]
                contig_path2 = path[:start_index]
                logger.info(f"hi conf node: {node}")
                # use the full contig score
                score = scores_dict[contig_id]
                # Path dont store the fisrt node
                # if score >= 0.5:
                path_dict[0].setdefault(contig_path1[0], []).append((contig_path1[1:], contig_id, score,contig_path2))
                logger.info(f"add contig path: {contig_path1[0]} {str(contig_path1[1:])}, score =  {score}, id={contig_id}")
                logger.info(f"It's pre contig is {contig_path2}")
            


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
    logger.info(f"merge {str(path1)} and {str(path2)}")
    right_idx-=1
    return path1[:left_idx] + path2[right_idx:]+path2[:right_idx] + path1[left_idx:]

def meet_criterion(path, G, SEQS,max_k, max_CV,bamfile):
    return get_wgtd_path_coverage_CV(path,G,SEQS,max_k_val=max_k) <= max_CV and is_good_cyc(path,G,bamfile)


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
