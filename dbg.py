"""

Author: Xisen Jin
Date: 2018.05.27

De Brujin Graph

"""

from gene import *
import collections

class Edge:
    def __init__(self, seq1, seq2, ori):
        # seq1 --> seq2
        self.seq1 = seq1
        self.seq2 = seq2
        self.ori = '' # '11','12','21','22'
        self.cov = 1

class Vertex:
    def __init__(self, seq):
        self.seq = seq
        self.cov = 1
        self.out_edge = {}

    def rev(self):
        return rev_complement(self.seq)

class DBG:
    """
    Data structure
    
    seq --O(1)--> vertex
    vertex + seq --O(1) --> edge
    
    """
    def __init__(self, k):
        self.k = k
        self.v = {}

    def add_vertex(self, kmer):
        kmer = canon(kmer) # store the smaller one
        if kmer in self.v:
            self.v[kmer].cov += 1
        else:
            self.v[kmer] = Vertex(kmer)

    def add_edge(self, kmer1, kmer2):
        """
        edge: kmer1 to kmer2
        """
        ckmer1 = canon(kmer1)
        ckmer2 = canon(kmer2)
        ori1 = '1' if ckmer1 == kmer1 else '2'
        ori2 = '1' if ckmer2 == kmer2 else '2'
        if ckmer2 not in self.v[ckmer1].out_edge:
            self.v[ckmer1].out_edge[ckmer2] = Edge(ckmer1, ckmer2, ori1 + ori2)
        else:
            self.v[ckmer1].out_edge[ckmer2].cov += 1
        if ckmer1 not in self.v[ckmer2].out_edge:
            self.v[ckmer2].out_edge[ckmer1] = Edge(ckmer2, ckmer1, ori2 + ori1)
        else:
            self.v[ckmer2].out_edge[ckmer1].cov += 1

    def add_read(self, read):
        if len(read) <= self.k:
            return
        kmer1, kmer2 = None, None
        for i in range(len(read)-self.k):
            kmer1 = read[i:i+self.k]
            kmer2 = read[i+1:i+self.k+1]
            self.add_vertex(kmer1)
            self.add_edge(kmer1, kmer2)
        self.add_vertex(kmer2)
        self.add_edge(kmer1, kmer2)

    def graph_simplification(self):
        """
        simplify the graph by merging nodes
        implemented with BFS
        """

        visit = {} # key is seq, not vertex

        def simplify_path(v, in_ori):
            """
            simplify a chain that starts from the v
            return a new vertex and final visited node
            """
            seqs = []
            final_edge = None
            while len(v.out_edge) <= 2:
                if not seqs:
                    seqs.extend(list(compif(v.seq, in_ori)))
                else:
                    seqs.append(compif(v.seq, in_ori)[0])
                visit[v] = True

                next_edge = None
                for edge in v.out_edge.values():
                    if in_ori == edge.ori[0] and edge.seq2 not in visit: # compatible
                        next_edge = edge
                        break
                if next_edge is None:
                    v = None
                    final_edge = None
                    break

                v = self.v[next_edge.seq2]
                final_edge = next_edge.seq2
                in_ori = edge.ori[1]
            
            return v, ''.join(seqs), final_edge

        def merge(s,e,new_seq, first_edge_seq,final_edge_seq):
            self.add_vertex(new_seq)
            s.out_edge.pop(first_edge_seq)
            self.add_edge(s.seq, new_seq)
            if e is not None:
                e.out_edge.pop(rev_complement(final_edge_seq))
                self.add_edge(new_seq, e.seq)

        while True:
            # initialize
            q = collections.deque()
            first_item = None
            for k in self.v:
                if len(self.v[k].out_edge) != 2 and not visit[k]:
                    first_item = self.v[k]
                    break
            if first_item is None:
                break
            q.append(first_item)

            while len(q):
                v = q.popleft()
                visit[v.seq] = True
                for edge in v.out_edge:
                    ori = edge.ori[1]
                    next_v, seq, final_edge_seq = simplify_path(v, ori)
                    merge(v, next_v, seq, v.seq, final_edge_seq)
                    if next_v is not None:
                        q.append(next_v)
        
    def output_edges():


        









