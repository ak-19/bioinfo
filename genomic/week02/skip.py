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