import os
import logging
from collections import defaultdict, deque
import math
from Bio import SeqIO
import heapq


#########################################
##########   Fasta handling   ##########
#########################################
def read_fasta_gen(path: str):
    """
    summary:
        Read a FASTA file and yield (record_id, sequence) pairs as a generator.

    complexity:
        Time: O(L), where L = total characters in FASTA.
        Memory: O(max_record_length); only one record stored at a time.

    description:
        1. Validate file extension (.fasta, .fa, .fna, .fas).
        2. Open FASTA file using Bio.SeqIO.
        3. For every record:
            - extract ID
            - convert sequence to uppercase
            - yield the pair (id, seq)
        4. If file not found, log error and raise FileNotFoundError.

    args:
        - path (str): path to FASTA file.

    returns:
        generator yielding (id, sequence)
    """
    
    basename = os.path.basename(path)
    
    if not any(path.endswith(ext) for ext in (".fasta", ".fa", ".fas", ".fna")):
        logging.error(f'Unsuported file extension for: {path}')
        raise ValueError(f'Given file: {basename} appears not to be in proper format - .fasta or .fna.\n full path: {path}\n')

    logging.info(f'Reading FASTA file: {basename}, full paht: {path}')
    
    try:
        for record in SeqIO.parse(path, "fasta"):
            yield record.id, str(record.seq).upper()
    except FileNotFoundError:
            logging.error((f'Reading fasta: Given file: {basename} appears not found.\n full path: {path}\n'))
            raise FileNotFoundError((f'Reading fasta: Given file: {basename} appears not found.\n full path: {path}\n'))
    



#########################################
####### Rolling/Rabin-KArp hasher #######
#########################################

class RabinKarpHash:
    """
    summary:
        Precompute prefix hashes and powers for O(1) substring hashing.

    time and memory complexity:
        Time: O(L) to build prefix and power arrays.
        Memory: O(L) for arrays.

    description:
        1. Validate input string (non-empty).
        2. Convert sequence to uppercase.
        3. Build pow[i] = base^i mod M.
        4. Build pref[i] = hash of s[0:i].
        5. Store arrays for O(1) substring extraction.

    args:
        - s (str): DNA sequence.
        - base (int): rolling hash base. - any number bigger that alphabet size - hash("ACG") = A * base^2 + C * base^1 + G * base^0
        Every letter has its weight depending on position
        - mod (int): a large prime modulus (2^61−1).   Using modular arithmetic reduces collisions and keeps hash values bounded.
        
    returns:
        RabinKarpHash object containing prefix hashes and powers.
    """
    
    def __init__(self, s: str, base: int = 127, mod: int = 2**61-1):
        if not isinstance(s, str):
            logging.error(f'RabinKarpHash: input is not a string')
            raise TypeError("Input to RabinKarpHash must be a string")
        
        if len(s) == 0:
            logging.error("RabinKarpHash: empty sequence provided")
            raise ValueError('Given DNA sequence is empty. Cannot create a hash')
                
        self.s = s.upper()
        self.n = len(self.s)
        self.base = base 
        self.mod = mod
        self.pow = [1] * (self.n + 1)
        for i in range(1, self.n + 1):
            self.pow[i] = (self.pow[i-1] * self.base) % self.mod
        self.pref = [0] * (self.n + 1)
        for i, ch in enumerate(self.s, start=1):
            v = self._ntoint(ch)
            self.pref[i] = (self.pref[i-1] * self.base + v) % self.mod



    @staticmethod
    def _ntoint(c: str) -> int:
        # map canonical A,C,G,T -> ints 1..4, unknown -> 5
        if c == 'A': return 1
        if c == 'C': return 2
        if c == 'G': return 3
        if c == 'T': return 4
        return 5
        

    def subhash(self, l: int, r: int) -> int:
        """
        summary:
            Return hash for substring s[l:r] in O(1).
            pref[i] is hash(s[0:i]) - its hash of every prefixed part of sequence, where i is a position

        time and memory complexity:
            Time: O(1)
            Memory: O(1)

        description:
            1. Validate substring indices.
            2. Compute hash using:
            pref[r+1] - pref[l] * pow[r-l+1]  (mod M)
            3. Return computed integer hash.

        args:
            - l (int): left index - first character inclusively
            - r (int): right index (inclusive) - last character inclusively

        returns:
            int: hash of substring
        """
        if l < 0 or r >= self.n or l > r:
            raise IndexError("Invalid substring indices")
        # pref indices are 1-based
        res = (self.pref[r+1] - self.pref[l] * self.pow[r - l + 1]) % self.mod
        return res

    
#########################################
####### Minimizer/winnowing index #######
#########################################
def compute_minimizer_index(reference: str, K: int = 15, W: int = 30, max_occ: int = 1000)-> dict:
    """
    summary:
        Build a minimizer index mapping minimizer_hash → list of reference positions (start of k-mer).
        Uses rolling hash - next hash is computed on the base of previous. Saves complexity
        
    time and memory complexity:
        Time: O(L) amortized; O(L*W) worst-case. Worst case occurs if the deque grows to size W and all subhashes are strictly increasing, causing each element to stay in the deque for W iterations before removal.
        Memory: O(L) positions stored; bounded per key by max_occ.

    description:
        1. Validate K and W (K <= W, lengths sufficient).
        2. Build rolling hash for reference.
        3. Slide a window of W over k-mer hashes.
        4. Maintain deque of candidate minimizers.
        5. For each full window:
            - select smallest hash as minimizer
            - append its position to dictionary entry if count < max_occ
        6. Return dictionary of minimizer → positions list.

    args:
        - reference (str): reference DNA
        - K (int): k-mer size
        - W (int): window size
        - max_occ (int): max stored occurrences

    returns:
        dict: minimizer_hash → [positions]
        Example output:
        {11494769: [7],
        54275408: [10],
        35984216: [39], ...
        Shows where window has the lowest hash
        
        keys are minimizers,- the smallest hash of k-mer in each window of reference
        values - each list contains the position in the reference, where that minimizer occurs

    """

    # Build rolling hash
    ref = reference.upper()
    L = len(ref)
    
    # -------- Validate input --------
    if K <= 0 or W <= 0 or K > W or K > L or W > L:
        raise ValueError("Bad K/W or too short reference")

    try:
        rk = RabinKarpHash(reference)
    except Exception as e:
        logging.error(f"Computing minimizer index: Failed to initialize RabinKarpHash: {e}")
        raise

    window_kmers = W - K + 1 # number of kmers in single window
    deq = deque() # begining = minimal hash in the window
    minimizer_index = defaultdict(list)
    
    for i in range(L - K + 1):
        h = rk.subhash(i, i + K - 1)
        # drop outside
        window_start = i - (window_kmers - 1)
        # delte kmers from outside the window - O(1) time
        while deq and deq[0][1] <= window_start:
            deq.popleft()
        # maintain increasing hashes: smallest at left. Delete biggger hashes from the end of que.
        while deq and deq[-1][0] >= h:
            deq.pop()
        # add current kmer
        deq.append((h, i))
        # save minimizer when window is full
        if i >= window_kmers - 1:
            mh, pos = deq[0]
            if len(minimizer_index[mh]) < max_occ:
                minimizer_index[mh].append(pos)
    return dict(minimizer_index)


#########################################
####### MinHash (bottom sketches) #######
######### for estimating Jaccar #########
#########################################

def build_kmer_hashes(seq: str, K: int):
    """
    summary:
        Yield rolling-hash values for all k-mers in the sequence.
        For incerase speed, use hash of forward only instead of reverse-complement (2x faster)

    time and memory complexity:
        Time: O(n) for n-k+1 k-mers.
        Memory: O(1).

    description:
        1. Convert sequence to uppercase.
        2. Build RabinKarpHash for the sequence.
        3. For each position i from 0 to n-K:
            - compute subhash(i, i+K-1)
            - yield hash value.

    args:
        - seq (str): DNA sequence
        - K (int): k-mer length

    returns:
        generator producing hashed k-mers
    """
    
    seq = seq.upper()
    n = len(seq)
    rk = RabinKarpHash(seq)
    for i in range(0, n - K + 1):
        ksub = seq[i:i+K]
        # choose min of forward and reverse complement string
        # But for speed we use hash of forward only
        yield rk.subhash(i, i+K-1)

def minhash_sketch_from_kmer_hashes(kmer_hash_iter, sketch_size: int = 200) -> set:
    """
    summary:
        Keep the sketch_size smallest hashes (bottom sketch).

    time and memory complexity:
        Time: O(n log s), where n is the number of k-mers, and s is the sketch_size.
        Memory: O(s) for heap and set.

    description:
        1. Validate sketch_size > 0.
        2. Create a max-heap (store negative values).
        3. For each hash h:
            - if heap not full → push -h
            - else if h < largest in heap → replace heap top
        4. Convert negative values back to positive and return set.

    args:
        - kmer_hash_iter: iterator of kmer hashes - yielding one hashed kmer at a time.
        - sketch_size (int): number of smallest hashes to keep. Representative probe of kmers with size sketch_size

    returns:
        set of int: A set containing the bottom sketch hashes.
    """
    
    if sketch_size <= 0:
        raise ValueError("sketch_size must be positive")
    # maintain a max-heap of size <= sketch_size (store negative to make max-heap)
    heap = []
    for h in kmer_hash_iter:
        if len(heap) < sketch_size:
            heapq.heappush(heap, -h) #add hashes untill heap is full
        else:
            if h < -heap[0]: # (meaning: new hash is smaller than the largest element in the current sketch).
                heapq.heapreplace(heap, -h)
    # return set of positive hashes
    return set([-x for x in heap])

def jaccard_from_sketches(sk1: set, sk2: set) -> float:
    """
    summary:
        Compute estimated Jaccard similarity using sketch intersection.

    time and memory complexity:
        Time: O(s) where s is sketch size (reminder: sketch is the probe of the kmers)
        Memory: O(s).

    description:
        1. If both sets empty, return 1.0.
        2. Compute intersection size.
        3. Compute union size.
        4. Return intersection / union.

    args:
        - sk1 (set): sketch of sequence A
        - sk2 (set): sketch of sequence B

    returns:
        float: Jaccard estimate
    """
    
    if not sk1 and not sk2:
        return 1.0
    inter = len(sk1.intersection(sk2))
    union = len(sk1) + len(sk2) - inter
    return inter / union if union > 0 else 0.0

#########################################
###### Reference window sketches #######
#########################################


def build_reference_window_sketches(reference: str, K: int, window_size: int, sketch_size: int = 200, stride: int = 100):
    """
    summary:
        Slide windows across reference and compute MinHash sketch for each window.

    time and memory complexity:
        Time: O((L / stride) * (window_size + sketch_size log sketch_size))
        Memory: O((L / stride) * sketch_size)

    description:
        1. For every window start (step = stride):
            - extract reference[start:end]
            - build k-mer hashes
            - compute bottom sketch
            - store (start, end, sketch)
        2. Ensure final window covers reference tail.
        3. Return list of window sketches.

    args:
        - reference (str): reference sequence
        - K (int): k-mer length
        - window_size (int): sliding window size
        - sketch_size (int): sketch size
        - stride (int): window step. Distance between another windows begginigs in reference

    returns:
        list of (start, end, sketch_set)
    """
    
    L = len(reference)
    windows = []
    
    for start in range(0, max(1, L - window_size + 1), stride):
        end = min(L, start + window_size)
        seq = reference[start:end]
        kiter = build_kmer_hashes(seq, K)
        sk = minhash_sketch_from_kmer_hashes(kiter, sketch_size)
        windows.append((start, end, sk))
        
    # ensure last window covers tail (f.e. if window is 200 length, ref has 1050 and strid = 100), last window is 900-1000
    if windows and windows[-1][1] < L:
        start = max(0, L - window_size)
        end = L
        seq = reference[start:end]
        sk = minhash_sketch_from_kmer_hashes(build_kmer_hashes(seq, K), sketch_size)
        windows.append((start, end, sk))
    return windows




#########################################
######  Seeding using minimizers  #######
#########################################

def seeding_from_minimizers(read: str, minimizer_index: dict, K: int, W: int, min_seeds: int = 2):
    """
    summary:
        Compute minimizers in read and map them to reference using minimizer index.

    time and memory complexity:
        Time: O(L_read + #minimizers * avg_occ)
        Memory: O(L_read) (for k-mer hash list); O(W) additional.

    description:
        1. Validate read length >= K and W.
        2. Compute k-mer hashes using RabinKarpHash.
        3. Slide window W to extract read minimizers.
        4. For each minimizer hash:
            - looking for minimizer_index[mh]
            - add (ref_pos, read_pos) to seed_hits
        5. Count distinct minimizers.
        6. If < min_seeds → unmapped.
        7. Return (seed_hits, unmapped)

    args:
        read (str): your query
        minimizer_index (dict): output of function: compute_minimizer_index().
        K (int): kmer length
        W (int): window size
        min_seeds (int, optional): minimum number of seeds to classify. Defaults to 2.

    returns:
        list: (seed_hits(list), unmapped(bool)) -> seed_hits: list of (ref_pos, read_pos)
    """
    
    read = read.upper()
    L = len(read)
    if L < K or L < W:
        return [], True
    # make hashes from read
    rk = RabinKarpHash(read)
    # split hashed read to kmers
    kmer_hashes = [rk.subhash(i, i+K-1) for i in range(L - K + 1)]
    window_kmers = W - K + 1
    
    deq = deque()
    read_minimizers = []
    for i in range(len(kmer_hashes)):
        while deq and deq[0] < i - window_kmers + 1:
            deq.popleft()
        while deq and kmer_hashes[deq[-1]] >= kmer_hashes[i]:
            deq.pop()
        deq.append(i)
        if i >= window_kmers - 1:
            mp = deq[0]
            mh = kmer_hashes[mp]
            read_minimizers.append((mh, mp))
            
    # map to reference positions
    seed_hits = []
    for mh, read_pos in read_minimizers:
        for ref_pos in minimizer_index.get(mh, []):
            seed_hits.append((ref_pos, read_pos))
    distinct = len(set([m for m, _ in read_minimizers]))
    unmapped = distinct < min_seeds #True or False
    return seed_hits, unmapped




#########################################
###############  Chaining  ##############
#########################################

def chain_seeds(seed_hits: list, cluster_gap: int = 200, min_cluster_seeds: int = 4):
    """
    summary:
        Group seed hits into clusters likely representing same alignment.

    time and memory complexity:
        Time: O(S log S), S = number of seed hits becouse standard python sort is S log s
        Memory: O(S)

    description:
        1. If no seeds, return empty list.
        2. Compute diagonal d = ref_pos - read_pos. For seeds coming from the same true alignment, ref_pos increases ~1 per base and read_pos increases ~1 per base, so their difference (diag) remains ~constant.
        3. Group seeds by diagonal buckets. Buckets are groups of seeds with similiar diag values.
        4. For each bucket:
            - sort seeds
            - greedily merge seeds into dense clusters. Greedy means: go through sorted list and extend current cluster if gap and diag conditions are done.
            - compute cluster stats (diag, start_ref, end_ref, score)
        5. Sort clusters by score.
        6. Return list of clusters.

    args:
        - seed_hits (list)
        - cluster_gap (int)
        - min_cluster_seeds (int)

    returns:
        list of cluster dicts:
        [
            {
                'diag': 990.0,            # average diag for first cluster
                'start_ref': 1000,
                'end_ref': 1051,         # end exclusive
                'score': 3,              # liczba seedów w klastrze
                'seeds': [(1000,10),(1010,20),(1050,60)]
            },
    """
    
    if not seed_hits:
        return []
    # group by coarse diagonal bucket
    buckets = defaultdict(list)
    for r, p in seed_hits:
        diag = r - p
        key = round(diag / max(1, cluster_gap))  # coarse bucket
        buckets[key].append((diag, r, p))
    clusters = []
    for key, seeds in buckets.items():
        seeds.sort(key=lambda x: (x[0], x[1]))  # by diag, ref_pos
        # greedy scan creating dense clusters based on ref_gap <= cluster_gap and diag closeness
        curr = [seeds[0]]
        for nxt in seeds[1:]:
            dprev, rprev, _ = curr[-1]
            dcur, rcur, _ = nxt
            if abs(rcur - rprev) <= cluster_gap and abs(dcur - dprev) <= 3 * cluster_gap:
                curr.append(nxt)
            else:
                if len(curr) >= min_cluster_seeds:
                    diags = [x[0] for x in curr]
                    refs = [x[1] for x in curr]
                    clusters.append({
                        'diag': sum(diags) / len(diags),
                        'start_ref': min(refs),
                        'end_ref': max(refs) + 1,  # end exclusive
                        'score': len(curr),
                        'seeds': [(x[1], x[2]) for x in curr]
                    })
                curr = [nxt]
        # final
        if len(curr) >= min_cluster_seeds:
            diags = [x[0] for x in curr]
            refs = [x[1] for x in curr]
            clusters.append({
                'diag': sum(diags) / len(diags),
                'start_ref': min(refs),
                'end_ref': max(refs) + 1,
                'score': len(curr),
                'seeds': [(x[1], x[2]) for x in curr]
            })
    # sort by score desc
    clusters.sort(key=lambda x: x['score'], reverse=True)
    return clusters


#########################################
########   Limiting candidates  #########
#########################################


def limit_candidates(candidates: list, max_cand: int = 5):
    """
    keeps only top-n clusters to increase algorithm speed

    Args:
        candidates (list): list of candidates - output of chain_seeds function
        max_cand (int, optional): Maximal number of candidats Defaults to 5.

    Returns:
        _type_: sorted by score list of candidates to be alignment place
    """
    
    """
    summary:
        Keep only the top-scoring alignment clusters.

    time and memory complexity:
        Time: O(C log C), C = number of clusters.
        Memory: O(C)

    description:
        1. If no candidates, return empty list.
        2. Sort candidates by score descending.
        3. Return top max_cand clusters.

    args:
        - candidates (list)
        - max_cand (int)

    returns:
        list of clusters
    """
    
    if not candidates:
        return []
    return sorted(candidates, key=lambda x: x['score'], reverse=True)[:max_cand]





#########################################
####  Computing banded edit distance ####
#########################################

def banded_edit_distance(read: str, ref: str, diag_center: int, max_error_ratio: float = 0.09):
    """
    summary:
        Compute banded edit distance around estimated diagonal.
        Instead of store whole DP in memory, we store just PREV and CURR rows. 

    time and memory complexity:
        Time: O(m * kmax), m = read length (=rows in DP)
        Memory: O(kmax)

    description:
        1. Convert read/ref to uppercase.
        2. Compute m (len read), n (len ref) and kmax = ceil(max_error_ratio * m). kmax is the maximum allowed deviation from the diagonal, proportional to allowed edit distance.
        3. Compute band width B = 2*kmax + 1. The band is a diagonal strip around the expected alignment path. Only DP cells with |(j - i) - diag_center| <= kmax are computed.
        4. Initialize DP row 0 with allowed columns.
        5. For each row i (1..m):
            - compute allowed j range
            - compute DP values (sub, del, ins)
            - track best alignment score and ref start
        6. Return best distance and start position.

    Args:
        read (str): your DNA query
        ref (str): reference DNA sequence
        diag_center (int): midpoint of possible place where read has to align to reference. Start point of searching perfect alignment
        max_error_ratio (float, optional): % indicator how much read can differ to alignment. Defaults to 0.15.

    Returns:
        _type_: If alignment is impossible returns (inf, None). Otherwise returns (best_dist, best_ref_start), where best_dist is minimal number of editions founded while DP. Best_ref_start is position in ref where the best alignemtn begins
    """
    
    read, ref = read.upper(), ref.upper()
    m, n = len(read), len(ref)
    if m == 0:
        return 0, 0
    kmax = math.ceil(max_error_ratio * m)
    B = 2 * kmax + 1
    INF = 10**9

    # i-th row is read[0:i] characters
    # j-th column is ref[0:j] characters
    prev = {}
    # initialize row i = 0 (empty prefix of read aligned to prefixes of ref)
    i = 0
    # allowed j range for i=0. If diag center position is faulty in length more tthan kmax, it wont find proper alignemnt
    jmin0 = max(0, i + diag_center - kmax)
    jmax0 = min(n, i + diag_center + kmax)
    for j in range(jmin0, jmax0 + 1):
        prev[j] = j  # cost of inserting j reference chars (i=0 -> need j insertions)

    best_dist = INF
    best_ref_start = None

    for i in range(1, m + 1):
        curr = {}
        jmin = max(0, i + diag_center - kmax)
        jmax = min(n, i + diag_center + kmax)
        if jmin > jmax:
            # band ran out of reference range
            break
        for j in range(jmin, jmax + 1):
            # substitution / match from (i-1, j-1)
            sub = INF
            if (j-1) in prev:
                cost = 0 if read[i-1] == ref[j-1] else 1
                sub = prev[j-1] + cost
            # deletion (from read) -> from (i-1, j)
            delete = prev.get(j, INF) + 1
            # insertion (into read) -> from (i, j-1)
            insert = curr.get(j-1, INF) + 1
            val = min(sub, delete, insert)
            curr[j] = val
        # track best in this row with j corresponding to end of ref alignment position
        row_min_j, row_min_val = min(curr.items(), key=lambda x: x[1])
        if row_min_val < best_dist:
            best_dist = row_min_val
            # compute where read[0] aligns: if at row i we are at ref position j, then alignment started at:
            # start = j - i
            best_ref_start = row_min_j - i
            if best_ref_start < 0:
                best_ref_start = 0
        prev = curr


    if best_dist >= INF:
        return INF, None
    return best_dist, best_ref_start 


#########################################
######### Verifying candidates ##########
#########################################
def verify_candidate(read: str, reference: str, candidate: dict, max_error_ratio: float = 0.20, margin: int = 50):
    """
    summary:
        Extract reference window around candidate and verify with banded DP.

    time and memory complexity:
        Time: O(m * kmax)
        Memory: O(kmax)

    description:
        1. Compute read length m and cluster diagonal.
        2. Estimate alignment start from cluster mid.
        3. Expand region by margin and kmax.
        4. Extract reference fragment.
        5. Run banded_edit_distance on fragment.
        6. Convert relative start to absolute start.
        7. Return {dist, start, end, score} dict.

    args:
        read (str): DNA seq
        reference (str): reference DNA seq
        candidate (dict): candidate to perfect alignment. 
        max_error_ratio (float, optional): Defaults to 0.20.
        margin (int, optional): Maximal bias to left and right  positions from diag. Defaults to 50.


    returns:
        dict: mapping result summary -dict: {'dist': int(dist), 'start': int(abs_start), 'end': int(abs_end), 'score': candidate.get('score', 0)}
        Score may represent seed count (supporting evidence for alignment).

    """
    
    m = len(read)
    # center diag
    diag = int(round(candidate['diag'])) #This is approximate ref_pos - read_pos for seeds.
    # estimate start of alignment in reference by centering read on cluster
    cluster_mid = (candidate['start_ref'] + candidate['end_ref']) // 2
    estimated_start = cluster_mid - (m // 2)
    
    # expand by margin and by kmax
    kmax = math.ceil(max_error_ratio * m)
    start = max(0, estimated_start - (kmax + margin))
    end = min(len(reference), estimated_start + m + kmax + margin)
    ref_fragment = reference[start:end]
    
    dist, rel_start = banded_edit_distance(read, ref_fragment, diag_center=diag - start, max_error_ratio=max_error_ratio)
    if rel_start is None: #relative start
        return {'dist': float('inf'), 'start': None, 'end': None, 'score': candidate.get('score', 0)}
    abs_start = start + rel_start
    abs_end = abs_start + m
    
    # Ensure thta computed end position doesnt go past the end of the reference
    abs_end = min(abs_end, len(reference))
    
    return {'dist': int(dist), 'start': int(abs_start), 'end': int(abs_end), 'score': candidate.get('score', 0)}


#########################################
##########  Ranking and output ##########
#########################################
def rank_and_output(read_id: str, read: str, verified_candidates: list, max_error_ratio=0.20):    
    """
    summary:
        Rank verified candidates and choose the best alignment.

    time and memory complexity:
        Time: O(C), C = number of verified candidates.
        Memory: O(1)

    description:
        1. Compute acceptable distance kaccept = floor(max_error_ratio * read_len).
        2. Filter candidates with dist <= kaccept and not None.
        3. If no valid candidate → return None.
        4. Rank by:
            - smallest distance
            - highest score
            - minimal length mismatch This penalizes candidates whose estimated aligned length differs from read length.
        5. Return [start, end] for best result.

    args:
        read_id (str): your query ID
        read (str): your DNA query
        verified_candidates (list): list of verified candidates. output of verify_candidate function
        max_error_ratio (float, optional): . Defaults to 0.20.


    returns:
        list [start, end] or None
    """
    
    m = len(read)
    if not verified_candidates:
        return None
    kaccept = math.floor(max_error_ratio * m)
    filtered = [c for c in verified_candidates if c['dist'] <= kaccept and c['start'] is not None]
    if not filtered:
        return None
    
    def rank_key(c):
        return (c['dist'], -c['score'], abs((c['end'] - c['start']) - m))
    best = min(filtered, key=rank_key)
    return [best['start'], best['end']]