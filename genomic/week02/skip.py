import bisect

class Index:
    def __init__(self, t, k):
        self.k = k
        self.index = []
        for i in range(len(t) + 1 - k):
            self.index.append((t[i:i+k], i))
        self.index.sort()
    
    def query(self, p):
        kmer = p[:self.k]
        
        start = bisect.bisect_left(self.index, (kmer, -1))
        
        result = []
        
        for i in range(start, len(self.index)):
            if self.index[i][0] == kmer:
                result.append(i)
            else:
                break
        
        return result

index = Index('AGAGACGATA', 3)
print(index.query('AGA'))


def boyer_moore(p, p_bm, t):

i = 0
occurrences = []
while i < len(t) - len(p) + 1:
shift = 1
mismatched = False
for j in range(len(p)-1, -1, -1):
if p[j] != t[i+j]:
skip_bc = p_bm.bad_character_rule(j, t[i+j])
skip_gs = p_bm.good_suffix_rule(j)
shift = max(shift, skip_bc, skip_gs)
mismatched = True
break
if not mismatched:
occurrences.append(i)
skip_gs = p_bm.match_skip()
shift = max(shift, skip_gs)
i += shift
return occurrences