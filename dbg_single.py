"""

Author: Xisen Jin
Date: 2018.05.27

De Brujin Graph

"""

from gene import *
import collections
import argparse

class Edge:
    def __init__(self, seq1, seq2):
        # seq1 --> seq2
        self.seq1 = seq1
        self.seq2 = seq2
        self.cov = 1

class Vertex:
    def __init__(self, seq):
        self.seq = seq
        self.cov = 1
        self.out_edge = {}
        self.in_edge = {}

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
        if kmer in self.v:
            self.v[kmer].cov += 1
        else:
            self.v[kmer] = Vertex(kmer)

    def add_edge(self, kmer1, kmer2):
        """
        edge: kmer1 to kmer2
        """
        if kmer2 not in self.v[kmer1].out_edge:
            self.v[kmer1].out_edge[kmer2] = Edge(kmer1, kmer2)
        else:
            self.v[kmer1].out_edge[kmer2].cov += 1
        if kmer1 not in self.v[kmer2].in_edge:
            self.v[kmer2].in_edge[kmer1] = Edge(kmer2, kmer1)
        else:
            self.v[kmer2].in_edge[kmer1].cov += 1

    def add_read(self, read):
        if len(read) <= self.k:
            return
        kmer1, kmer2 = None, None
        for i in range(len(read)-self.k):
            kmer1 = read[i:i+self.k]
            kmer2 = read[i+1:i+self.k+1]
            self.add_vertex(kmer1)
            self.add_vertex(kmer2)
            self.add_edge(kmer1, kmer2)

            self.add_vertex(rev_complement(kmer1))
            self.add_vertex(rev_complement(kmer2))
            self.add_edge(rev_complement(kmer2), rev_complement(kmer1))


    def dfs_wrapper(self, start_node, visit, ignore_low_cov=False):
        path = list(start_node.seq[:-1])
        node_stack = collections.deque()
        node_stack.append(start_node)
        while len(node_stack):
            node = node_stack.pop()
            path.append(node.seq[-1])
            if node.seq not in visit:
                visit[node.seq] = 1
            else:
                visit[node.seq] += 1
            for edge in node.out_edge.values():
                if not ignore_low_cov:
                    if edge.seq2 not in visit or self.v[edge.seq2].cov > visit[edge.seq2]:
                        node_stack.append(self.v[edge.seq2])

        return path

    def get_contig_wrapper(self, start_node, visit):

        def is_self_ring(v):
            return v.seq in v.out_edge

        def select_valid_path(v):
            """
            select a valid path to traverse
            """
            if is_self_ring(v):
                return None
            cmpt_edge = list(v.out_edge.values())
            # un-ambiguous path
            if len(cmpt_edge) == 1:
                return cmpt_edge[0]
            
            elif len(cmpt_edge) >= 2:
                cmpt_edge = sorted(cmpt_edge, key=lambda x:x.cov)
                if len(v.in_edge) == 1:
                    return cmpt_edge[-1]
                       

        node_stack = collections.deque()
        node_stack.append(start_node)
        contigs = []
        a = 0
        # stack stores nodes whose degree != 2
        while len(node_stack):
            a += 1
            b = 0
            node = node_stack.pop()
            visit.add(node.seq)
            flg = True
            #print(len(self.v), len(visit), a)
            for edge in node.out_edge.values():
                # travel until degree > 2 or degree = 1
                prefix = ''
                path = []
                v = self.v[edge.seq2]
                while True:
                    b += 1
                    path.append(v.seq[-1])
                    if v.seq in visit:
                        break
                    visit.add(v.seq)
                    sel_edge = select_valid_path(v)
                    if sel_edge is None:
                        if len(v.out_edge) >= 2:
                            node_stack.append(v)
                        break
                    else:
                        v = self.v[sel_edge.seq2]
                contigs.append(prefix + ''.join(path))
                flg = False
        return contigs           


    def dfs_graph(self):
        visit = {}
        all_path = []
        c = 0

        for k in self.v:
            if k not in visit and len(self.v[k].out_edge)!=2:
                start_node = self.v[k]
                path = self.dfs_wrapper(start_node, visit)
                all_path.append(''.join(path))
        all_path.sort(key=lambda x:-len(x))
        return all_path

    def get_contig_graph(self):
        """
        dfs and store any unambiguous contig
        """
        visit = set()
        all_path = []
        c = 0
        for v in self.v:
            if v not in visit and len(self.v[v].out_edge) != 2:
                start_node = self.v[v]
                paths = self.get_contig_wrapper(start_node, visit)
                all_path.extend(paths)
                c += 1
                print(c,len(self.v),len(visit))
        return all_path

def write_fa(f, lines):
    lines.sort(key=lambda x: len(x))
    for i,line in enumerate(lines):
            if len(line.strip()) > 50:
                f.write('>{} length {} xxx\n'.format(i*2+1, len(line.strip())))
                f.write(line.strip()+'\n')
                #f.write('>{} length {} xxx\n'.format(i*2+1, len(line.strip())))
                #f.write(rev_complement(line.strip())+'\n')


def run(kk):
    g = DBG(kk)
    reads = read_fasta('data/data1/short_1.fasta')
    reads += read_fasta('data/data1/short_2.fasta')
    #reads = ['TTCTACTATCGCTGTGGGATGGATCATAAA',
           #  'TTCTACTATCGCTGTGGGATGGATCATCCC',
    #         'AAATTCCCCCCCCCCCCCC']
    for i,read in enumerate(reads):
        g.add_read(read)
        print(i)
        #if i == 100:
        #    break
    contigs = g.get_contig_graph()
    #contigs = g.dfs_graph()
    f = open('result/data1_contig_%d.txt' % kk,'w')
    write_fa(f, contigs)

if __name__ == '__main__':
    run(17)
    #for kk in range(51,79,2):
    #    run(kk)