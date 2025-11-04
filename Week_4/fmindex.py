# >>===================================================<<
# ||+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+||
# ||--------------------|B|O|X|-----------------------|||
# |||B|u|r|r|o|w|s|-|W|h|e|e|l|e|r|-|T|r|a|n|s|f|o|r|m|||
# ||+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+||
# >>===================================================<<

# Template by Leo Ackermann (2025)
# Code by ____

from ks import simple_kark_sort



# >>==========================================================================<<
# >>                              FM-index class                              <<
# >>==========================================================================<<

class FMindex:

    # >>===[ Class constructor ]================================================

    def __init__(self, seq, verbose=False):
        self.seq = seq
        self.verbose = verbose

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



    # >>===[ Attribute initialisation functions ]===============================

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
            raise ValueError("BWT non calculée.")
        
        #Compter les fréquences totales de chaque caractère
        total_counts = {}
        for char in self.bwt:
            if char in total_counts:
                total_counts[char] += 1
            else:
                total_counts[char] = 1
            
        #Calculer fm_count
        sorted_chars = sorted(total_counts.keys())
        self.fm_count = {}
        cumulative_count = 0
        for char in sorted_chars:
            self.fm_count[char] = cumulative_count
            cumulative_count += total_counts[char]

        #Construire la table des rangs
        self.fm_rank = {}
        for ch in sorted_chars:
            self.fm_rank[ch] = [0] * (len(self.bwt) + 1)

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


    # >>===[ Question 2.2 functions ]===========================================

    def get_string(self):
        if self.bwt is None:
            raise ValueError("BWT non calculée.")
        if self.fm_count is None or self.fm_rank is None:
            raise ValueError("FM-index non initialisé.")
        
        n = len(self.bwt)
        result = [""] * n
        idx = 0
        for i in range(n):
            if self.bwt[i] == '$':
                idx = i
                break

        for i in range(n - 1, -1, -1):
            c = self.bwt[idx]
            result[i] = c
            idx = self.fm_count[c] + self.fm_ranks[c][idx]
        
        return ''.join(result).rstrip('$')
    


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
    

    # == 2.5 : Approximate pattern matching (≤ e erreurs) ==
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

# >>===[ Function main ]===========================================
if __name__ == "__main__":

    fm = FMindex("banana", verbose=True)

    print("\n" + "=" * 10 + "Reconstruction du string" + "=" * 10 + "\n")
    print("Reconstructed string (naive):", fm.get_string__naive())
    print("Reconstruction string (rapide):", fm.get_string())

    print("\n" + "=" * 10 + "Compression de la BWT" + "=" * 10 + "\n")

    encoded = fm.compress_rle_bwt()
    decoded = fm.decompress_rle_bwt(encoded)

    print(f"RLE encodé : : {encoded}")
    print(f"\nRLE décodé :  {decoded}")
    print("\nMTF Encodé :", fm.compress_mtf_bwt())
    
    print("\n" + "=" * 10 + "Recherche de motifs exact" + "=" * 10 + "\n")

    motifs = ["ana", "na", "ban", "a$", "z"]
    for m in motifs:
        print(f" Motif '{m}':")
        print(" - Présent :  ", fm.membership(m))
        print(" - Nb occur :  ", fm.count(m))
        print(" - Positions : ", fm.locate(m))

    print("\n" + "=" * 10 + "Recherche de motifs exact" + "=" * 10 + "\n")

    motifs = ["ana", "ann", "ban", "bana", "nana"]
    for m in motifs:
        print(f"Motif '{m}' (error e=1):")
        print(" - Position approx :", fm.locate_approx(m, e=1))