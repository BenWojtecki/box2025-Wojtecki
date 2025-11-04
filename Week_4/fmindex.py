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
            raise ValueError("BWT has not been computed yet.")
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
                raise ValueError("BWT has not been computed yet.")
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
            raise ValueError("BWT has not been computed yet.")
        
        #Compter les fréquences totales de chaque caractère
        total_counts = {}
        for char in self.bwt:
            if char in total_counts:
                total_counts[char] += 1
            else:
                total_counts[char] = 1
            
        #Calculer nombre total de caractères plus petits que c
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

# >>===[ Function main ]===========================================
if __name__ == "__main__":
    fm = FMindex("banana", verbose=True)
    print("Reconstructed string (naive):", fm.get_string__naive())
    rle_encoded = fm.compress_rle_bwt()
    rle_decoded = fm.decompress_rle_bwt(rle_encoded)
    print("\nRLE Encodé :", rle_encoded)
    print("RLE Décodé :", rle_decoded)
    print("\nMTF Encodé :", fm.compress_mtf_bwt())
    print("\nInversion rapide (LF-mapping):", fm.get_string())
