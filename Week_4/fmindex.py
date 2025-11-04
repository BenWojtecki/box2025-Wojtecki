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

        if self.verbose:
            print("Suffix array :", self.sa)
            print("BWT          :", self.bwt)



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

    # >>===[ Pattern matching functions ]=======================================




# >>===[ Function main ]===========================================
if __name__ == "__main__":
    fm = FMindex("banana", verbose=True)
    print("Reconstructed string (naive):", fm.get_string__naive())
    rle_encoded = fm.compress_rle_bwt()
    rle_decoded = fm.decompress_rle_bwt(rle_encoded)
    print("\nRLE Encodé :", rle_encoded)
    print("RLE Décodé :", rle_decoded)
    print("\nMTF Encodé :", fm.compress_mtf_bwt())
