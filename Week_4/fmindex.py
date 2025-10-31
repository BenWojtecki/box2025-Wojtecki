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
        self.sa = None
        self.bwt = None 
        self.fm_count = None
        self.fm_rank = None
        self.fm_ranks = None
        self.next_smallest_letter = None



    # >>===[ Attribute initialisation functions ]===============================

    def __compute_sa(self, seq):
        self.sa = simple_kark_sort(seq + '$')

    def __compute_bwt_from_sa(self, seq):
        self.bwt = ''.join(map(lambda x: seq[x-1] if x > 0 else '$', self.sa))

    def get_string__naive(self):
        first_cols = ["" for _ in range(len(self.bwt))]
        for _ in range(len(self.bwt)):
            first_cols = sorted(map(lambda x: x[0] + x[1], zip(self.bwt, first_cols)))

        return first_cols[0][1:]




    # >>===[ Pattern matching functions ]=======================================



