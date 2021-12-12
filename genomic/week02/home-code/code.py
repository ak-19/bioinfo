import bisect

from bm_preproc import BoyerMoore

def naive_with_counts(pattern, text):
    P, T = len(pattern), len(text)
    occurence = []
    
    character_comparisons = 0 
    alignments_tried = 0
    
    for i in range(T - P + 1):
        alignments_tried += 1
        j = 0 
        while j < P:
            character_comparisons += 1
            if pattern[j] != text[j + i]: 
                break
            j += 1
            
        if j == P: occurence.append(i)
    
    return occurence, alignments_tried, character_comparisons

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

genome = readGenome('019-chr1.GRCh38.excerpt.fasta')


def boyer_moore_with_counts(p, p_bm, t):
    i = 0
    occurrences = []
    alignments_tried, character_comparisons = 0, 0
    while i < len(t) - len(p) + 1:
        alignments_tried += 1
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            character_comparisons += 1
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
        
    return occurrences, alignments_tried, character_comparisons


class Index(object):
    """ Holds a substring index for a text T """

    def __init__(self, t, k):
        """ Create index from all substrings of t of length k """
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    def query(self, p):
        """ Return index hits for first k-mer of p """
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

def app_match(p, t, n):
    segment_len = round(len(t) / (n+1))
    matches = set()
    for i in range(n + 1):
        start = i * segment_len
        end = min((i + 1) * segment_len, len(p))
        bm = BoyerMoore(p[start:end], alphabet='ACGT')
        curr_matches, a, b = boyer_moore_with_counts(p[start:end], bm, t)
        for m in curr_matches:
            if m < start or m-start+len(p) > len(t): continue
                
            misses = 0
            for j in range(0, start):
                if p[j] != p[m-start+j]:
                    misses += 1
                    if misses > n:
                        break
            
            for j in range(end, len(p)):
                if p[j] != p[m-start+j]:
                    misses += 1
                    if misses > n:
                        break                    
            
            if misses <= n:
                matches.add(m-start)
        
        return list(matches), a, b

# print(app_match('GGCGCGGTGGCTCACGCCTGTAAT', genome, 2))
# print(naive_with_counts('GGCGCGGTGGCTCACGCCTGTAAT', genome)   )

class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

def query_subseq(p, t, ind):
    n = 2
    all_matches = set()
    num_index_hits = 0
    P = len(p)
    T = len(t)
    for i in range(P - ind.k + 1):
        cut = p[i:]
        matches = ind.query(cut)
        # Extend matching segments to see if whole p matches
        for m in matches:
            num_index_hits += 1
            mismatches = 0
            for j in range(P):
                if not p[j] == t[m+j]:
                    mismatches += 1
                    if mismatches > n:
                        break

            if mismatches <= n:
                all_matches.add(m)
    return list(all_matches), num_index_hits

t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
p = 'to-morrow and to-morrow '
subseq_ind = SubseqIndex(t, 8, 3)

occurrences, num_index_hits = query_subseq(p, t, subseq_ind)

print(occurrences, num_index_hits)