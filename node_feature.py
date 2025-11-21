import numpy as np
import itertools
from collections import defaultdict

def dna_bases(k):
    """Generate all possible DNA k-mers in sorted order."""
    return [''.join(p) for p in itertools.product('ACGT', repeat=k)]

# 全局互补表（只创建一次）
_COMPLEMENT_TABLE = str.maketrans('ACGT', 'TGCA')  # 注意顺序对应！

def reverse_complement(kmer: str) -> str:
    """Return the reverse complement of a DNA k-mer (fast version)."""
    return kmer.translate(_COMPLEMENT_TABLE)[::-1]

def canonical_kmer(kmer: str) -> str:
    """Return the canonical k-mer (lexicographically smaller of kmer and its RC)."""
    rc = reverse_complement(kmer)
    return kmer if kmer < rc else rc

def get_kmers(seq, k):
    """Generate all k-mers from a DNA sequence."""
    for i in range(len(seq) - k + 1):
        yield seq[i:i+k]

def contig_to_freq_vector(contig_seq, k=4):
    """
    Convert a contig DNA sequence to a normalized k-mer frequency vector 
    of length 136, where k-mers are merged with their reverse complements.
    
    Returns a fixed-length numpy vector (length = 136), ordered by canonical k-mer.
    
    Args:
        contig_seq (str): DNA sequence of the contig (e.g., "ATGCGAT...")
        k (int): k-mer length (default: 4)
    
    Returns:
        np.ndarray: shape (136,) - frequency vector of canonical 4-mers
    """
    if len(contig_seq) < k:
        return np.zeros(136, dtype=np.float32)
    
    # Step 1: Count k-mers, grouped by canonical form
    canon_counts = defaultdict(int)
    seq = contig_seq.upper()
    
    for kmer in get_kmers(seq, k):
        # Only non-IUPAC bases
        if any(b not in 'ACGT' for b in kmer):
            continue
        canon = canonical_kmer(kmer)
        canon_counts[canon] += 1

    # Step 2: Total number of valid k-mers
    total_kmers = len(seq) - k + 1
    if total_kmers == 0:
        return np.zeros(136, dtype=np.float32)

    # Step 3: Get all canonical 4-mers (unique equivalence classes)
    all_kmers = dna_bases(k)
    canonical_forms = sorted(set(canonical_kmer(kmer) for kmer in all_kmers))
    
    # Verify dimension
    assert len(canonical_forms) == 136, f"Expected 136 canonical 4-mers, got {len(canonical_forms)}"

    # Step 4: Build 136-dimensional vector
    freq_vector = np.zeros(136, dtype=np.float32)
    for idx, canon in enumerate(canonical_forms):
        count = canon_counts.get(canon, 0)
        freq_vector[idx] = count / total_kmers  # normalized frequency

    return freq_vector  # shape: (136,)


def contig_set_to_freq_vector(contig_set: set, k=4):
    """
    Convert a set of contig DNA sequences into a single normalized k-mer frequency vector
    of length 136 (for k=4), using canonical k-mers (merged with reverse complements).

    Args:
        contig_set (set of str): Set of DNA sequences (e.g., {"ATGC...", "CGTA..."})
        k (int): k-mer length (default: 4)

    Returns:
        np.ndarray: shape (136,) - normalized frequency vector over canonical 4-mers
    """
    if not contig_set:
        return np.zeros(136, dtype=np.float32)

    # Step 1: Aggregate k-mer counts across all contigs
    canon_counts = defaultdict(int)
    total_kmers = 0

    for seq in contig_set:
        seq = seq.upper()
        if len(seq) < k:
            continue
        n_kmers_in_seq = len(seq) - k + 1
        valid_kmers_in_seq = 0

        for kmer in get_kmers(seq, k):
            if any(b not in 'ACGT' for b in kmer):
                continue
            canon = canonical_kmer(kmer)
            canon_counts[canon] += 1
            valid_kmers_in_seq += 1

        total_kmers += valid_kmers_in_seq  # only count valid ACGT kmers

    # Step 2: Handle edge case
    if total_kmers == 0:
        return np.zeros(136, dtype=np.float32)

    # Step 3: Generate canonical k-mer ordering (must match contig_to_freq_vector)
    all_kmers = dna_bases(k)
    canonical_forms = sorted(set(canonical_kmer(kmer) for kmer in all_kmers))
    assert len(canonical_forms) == 136, f"Expected 136 canonical 4-mers, got {len(canonical_forms)}"

    # Step 4: Build frequency vector
    freq_vector = np.zeros(136, dtype=np.float32)
    for idx, canon in enumerate(canonical_forms):
        count = canon_counts.get(canon, 0)
        freq_vector[idx] = count / total_kmers

    return freq_vector