import argparse
import gzip
from readfa import readfq
from ks import simple_kark_sort

class FMindex:

    # >>===[ Class constructor ]================================================

    def __init__(self, seq, verbose=False, k=20):
        self.seq = seq
        self.verbose = verbose
        self.k = k

        self.sa = None
        self.bwt = None 
        self.fm_count = None
        self.fm_rank = None
        self.fm_ranks = None
        self.next_smallest_letter = None

        self.__compute_sa(seq)
        self.__compute_bwt_from_sa(seq)
        self.__build_fm_index()

        if self.verbose:
            print("Suffix array :", self.sa)
            print("BWT          :", self.bwt)
            print("FM count     :", self.fm_count)
            print("FM rank (partial):")
            for k in self.fm_rank:
                print(f"   {k}: {self.fm_rank[k][:10]}")

    # >>===[ Question 1.1 functions ]===========================================
    def __compute_sa(self, seq):
        self.sa = simple_kark_sort(seq + '$')

    def __compute_bwt_from_sa(self, seq):
        self.bwt = ''.join(map(lambda x: seq[x-1] if x > 0 else '$', self.sa))

    # >>===[ Question 1.2 functions ]===========================================
    def get_string__naive(self):
        first_cols = ["" for _ in range(len(self.bwt))]
        for _ in range(len(self.bwt)):
            first_cols = sorted(map(lambda x: x[0] + x[1], zip(self.bwt, first_cols)))

        return first_cols[0][1:]

    # >>===[ Question 1.3 functions ]===========================================
    def compress_rle_bwt(self):
        if self.bwt is None:
            raise ValueError("BWT non calculé.")
        text = self.bwt

        if len(text) == 0:
            return ""
        
        encoded = ""
        count = 1
        for i in range(1, len(text)):
            if text[i] == text[i - 1]:
                count += 1
            else:
                encoded += text[i - 1] + str(count)
                count = 1

        encoded += text[-1] + str(count)
        if self.verbose:
            print("RLE Compressed BWT:", encoded)
        return encoded
    
    def decompress_rle_bwt(self, encoded_str):
        
        if(len(encoded_str) == 0):
            return ""
        
        decoded = ""
        i = 0
        
        while i < len(encoded_str):
            char = encoded_str[i]
            i += 1
            count_str = ""

            while i < len(encoded_str) and encoded_str[i].isdigit():
                count_str += encoded_str[i]
                i += 1

            count = int(count_str)
            for _ in range(count):
                decoded += char

        if self.verbose:
            print("Decompressed BWT:", decoded)
        return decoded
    

    def compress_mtf_bwt(self, text=None):
        if text is None:
            if self.bwt is None:
                raise ValueError("BWT non calculé.")
            text = self.bwt 
     
        alphabet = sorted(set(text))
        output = []

        for char in text:
            #Cherche la position de char dans alphabet
            pos = 0
            while pos < len(alphabet) and alphabet[pos] != char:
                pos += 1
            output.append(pos)

            saved = alphabet[pos]
            j = pos
            while j > 0:
                alphabet[j] = alphabet[j - 1]
                j -= 1
            alphabet[0] = saved 

        if self.verbose:
            print("MTF Compressed BWT:", output)
        return output

    # >>===[ Question 2.1 functions ]=======================================
    def __build_fm_index(self):
        if self.bwt is None:
            raise ValueError("BWT non calculé.")
        
        #Compter les fréquences totales de chaque caractère
        total_counts = {}
        for char in self.bwt:
            total_counts[char] = total_counts.get(char, 0) + 1
            
        #Calculer fm_count
        sorted_chars = sorted(total_counts.keys())
        self.fm_count = {}
        cumulative_count = 0
        for char in sorted_chars:
            self.fm_count[char] = cumulative_count
            cumulative_count += total_counts[char]

        #Construire la table des rangs
        self.fm_rank = {ch: [0] * (len(self.bwt) + 1) for ch in sorted_chars}

        for i in range(1, len(self.bwt) + 1):
            char = self.bwt[i - 1]
            for ch in sorted_chars:
                self.fm_rank[ch][i] = self.fm_rank[ch][i - 1]
            self.fm_rank[char][i] += 1
    
        self.fm_ranks = self.fm_rank

        #Construire next_smallest_letter
        self.next_smallest_letter = {}
        prev = None

        for ch in sorted_chars:
            self.next_smallest_letter[ch] = prev
            prev = ch

    # >>>>===[ Question 2.3 functions ]===========================================
    def __backward_search(self, pattern):
        l = 0
        r = len(self.bwt)
        for i in range(len(pattern) - 1, -1, -1):
            c = pattern[i]
            if c not in self.fm_count:
                return (0, 0)
            l = self.fm_count[c] + self.fm_rank[c][l]
            r = self.fm_count[c] + self.fm_rank[c][r]
            if l >= r:
                return (0, 0)
        return (l, r)

    def membership(self, pattern):
        l, r = self.__backward_search(pattern)
        return l < r

    def count(self, pattern):
        l, r = self.__backward_search(pattern)
        return r - l

    def locate(self, pattern):
        l, r = self.__backward_search(pattern)
        positions = []
        for i in range(l, r):
            positions.append(self.sa[i])
        positions.sort()
        return positions
    
    def locate_approx(self, pattern, e=1):
        results = []
        self.__approx_search_recursive(pattern, len(pattern) - 1, 0, len(self.bwt), e, results)

        positions = []
        for (l, r) in results:
            for i in range(l, r):
                positions.append(self.sa[i])

        unique_positions = []
        for p in sorted(set(positions)):
            unique_positions.append(p)
        return unique_positions
    
    def __approx_search_recursive(self, pattern, i, l, r, e, results):
        if e < 0:
            return

        if i < 0:
            if l < r:
                results.append((l, r))
            return

        current_char = pattern[i]

        # Cas même caractère
        if current_char in self.fm_count:
            new_l = self.fm_count[current_char] + self.fm_rank[current_char][l]
            new_r = self.fm_count[current_char] + self.fm_rank[current_char][r]
            if new_l < new_r:
                self.__approx_search_recursive(pattern, i - 1, new_l, new_r, e, results)

        # Cas erreur =>(substitution)
        if e > 0:
            for alt in self.fm_count.keys():
                if alt == current_char:
                    continue
                alt_l = self.fm_count[alt] + self.fm_rank[alt][l]
                alt_r = self.fm_count[alt] + self.fm_rank[alt][r]
                if alt_l < alt_r:
                    self.__approx_search_recursive(pattern, i - 1, alt_l, alt_r, e - 1, results)

def edit_distance(seq1, seq2, max_dist=None):
    n, m = len(seq1), len(seq2)

    if (max_dist is not None and abs(n - m) > max_dist):
        return max_dist + 1

    prev = list(range(m + 1))
    curr = [0] * (m + 1)

    for i in range(1, n + 1):
        curr[0] = i
        row_min = curr[0]
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                cost = 0
            else:
                cost = 1

            ins = curr[j-1] + 1
            dele = prev[j] + 1
            sub = prev[j-1] + cost

            v = ins if ins < dele else dele
            v = v if v < sub else sub

            curr[j] = v
            if v < row_min:
                row_min = v

        if (max_dist is not None and row_min > max_dist):
            return max_dist + 1

        prev, curr = curr, prev

    return prev[m]

def edit_ops(seq1, seq2):
    n, m = len(seq1), len(seq2)

    dp = [[0] * (m + 1) for _ in range(n + 1)]
    bt = [[None] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = i
        bt[i][0] = 'D'

    for j in range(1, m + 1):
        dp[0][j] = j
        bt[0][j] = 'I'

    for i in range(1, n+1):
        for j in range(1, m+1):
            if seq1[i-1] == seq2[j-1]:
                cost = 0
            else:
                cost = 1
            best = min(
                (dp[i-1][j] + 1, 'D'),
                (dp[i][j-1] + 1, 'I'),
                (dp[i-1][j-1] + cost, '=' if cost == 0 else 'X')
            )
            dp[i][j], bt[i][j] = best

    ops = []
    i, j = n, m
    while i > 0 or j > 0:
        op = bt[i][j]
        ops.append(op)
        if op == 'D':
            i -= 1
        elif op == 'I':
            j -= 1
        else:
            i -= 1
            j -= 1

    return dp[n][m], "".join(reversed(ops))

def revcomp(seq):
    comp_dict = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
        'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'
    }

    comp_seq = []
    for base in seq:
        comp_seq.append(comp_dict.get(base, 'N'))
    return ''.join(comp_seq)[::-1]

def get_kmer_seeds(query, fm, k):
    candidates = set()
    n = len(query)
    if n < k:
        return []

    for i in range(n - k + 1):
        kmer = query[i:i+k]
        l, r = fm._FMindex__backward_search(kmer)
        if l < r:
            for j in range(l, r):
                pos = fm.sa[j] - i
                if pos >= 0:
                    candidates.add(pos)

    return list(candidates)

def get_pos(hit):
    return hit[0]

def try_align(query, ref, fm, strand, max_edit=3, cluster_window=20):
    seed_positions = get_kmer_seeds(query, fm, fm.k)
    candidates = set(seed_positions)

    if not candidates:
        candidates = set(fm.locate_approx(query, e=min(max_edit, 2)))
    if not candidates:
        candidates = set(fm.locate(query))

    hits = []
    for pos in sorted(candidates):
        if pos < 0 or pos >= len(ref):
            continue
        ref_sub = ref[pos : pos + len(query)]
        if len(ref_sub) < len(query):
            continue
        d, ops = edit_ops(query, ref_sub)
        if d <= max_edit:
            hits.append((pos, d, ops, strand))

    if not hits:
        return []

    hits.sort(key=get_pos)
    merged = []
    current = hits[0]
    for h in hits[1:]:
        if abs(h[0] - current[0]) <= cluster_window:
            if h[1] < current[1]:
                current = h
        else:
            merged.append(current)
            current = h
    merged.append(current)

    return merged

def align_read_with_fm(read, ref, fm, max_edit=3):
    forward_hits = try_align(read, ref, fm, '+', max_edit=max_edit)
    reverse_hits = try_align(revcomp(read), ref, fm, '-', max_edit=max_edit)

    all_hits = forward_hits + reverse_hits

    if not all_hits:
        return None

    best = min(all_hits, key=lambda x: x[1])
    return best

def main():
    parser = argparse.ArgumentParser(description="FM-index based aligner")
    parser.add_argument("-k", type = int, default = 20, help="k-mer size (default 20)")
    parser.add_argument("ref")
    parser.add_argument("reads")
    parser.add_argument("-e", type = int, default=2, help="max edit distance, (default 2)")
    args = parser.parse_args()

    with open(args.ref) as f:
        ref_dict = {name: seq for name, seq, _ in readfq(f)}

    if len(ref_dict) == 0:
        print("No reference sequences found")
        return

    fm_indexes = {}
    for name, seq in ref_dict.items():
        fm_indexes[name] = (seq, FMindex(seq, k=args.k))

    with open(args.reads) as rf:
        reads_iter = readfq(rf)

        print("{:40s}	{:15s}	{:10s}	{:6s}	{:8s}	{}".format("#read", "ref", "pos", "strand", "score", "ops"))

        for rname, rseq, _ in reads_iter:
            pos = None
            dist = None
            ops = None
            strand = None
            best = None

            if len(rname) > 40:
                short_name = rname[:38] + "…"
            else:
                short_name = rname

            for ref_name, (ref_seq, fm) in fm_indexes.items():
                res = align_read_with_fm(rseq, ref_seq, fm, max_edit=args.e)

                if res is None:
                    continue

                # res is (pos, dist, ops, strand)
                if best is None or res[1] < best[2]:
                    best = (ref_name, res[0], res[1], res[2], res[3])

            if best is None:
                print("{:40s}	{:15s}	{:10s}	{:6s}	{:8s}	-".format(short_name, "-", "-", "-", "-"))
            else:
                ref_name, pos, dist, ops, strand = best
                score = len(rseq) - dist
                print("{:40s}	{:15s}	{:<10d}	{:6s}	{:<8d}	{}".format(short_name, ref_name, pos, strand, score, ops))


if __name__ == "__main__":
    main()
