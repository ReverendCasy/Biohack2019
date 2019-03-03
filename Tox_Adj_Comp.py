from Bio import SeqIO
from Bio import AlignIO
import numpy as np
from collections import defaultdict
from graphviz import Graph
from collections import deque
from itertools import chain
import pandas as pd
import argparse
#import subprocess


class AAGraph():
    def __init__(self, alignment):
        self.alignment = open(alignment)
        self.al = AlignIO.read(self.alignment, 'fasta')
        self.adjMat = np.triu(np.full((len(self.al), len(self.al)), 1))
        np.fill_diagonal(self.adjMat, 0)
        self.distMat = np.zeros((len(self.al), len(self.al)))
        self.alignment.close()
        self.vertices = list()
        self.edges = defaultdict(list)

    def upload_vertices(self):
        for entry in self.al:
            self.vertices.append(entry.id)

    def matBuild(self):
        self.distMat = np.zeros((len(self.al), len(self.al)))
        for k in range(1, len(self.al)):
            for i in range(k):
                for j in range(0, len(self.al[i].seq)):
                    if self.al[k].seq[j] != self.al[i].seq[j]:
                        self.distMat[i][k] += 1
                        self.edges[(self.al[k].id, self.al[i].id)].\
                            append(f"{self.al[k].seq[j]}{str(j + 1)}{self.al[i].seq[j]}")
        self.distMat = self.distMat

    def TR(self):
        for k in range(1, len(self.distMat)):
            for i in range(k):
                for j in range(k + 1, len(self.distMat)):
                    viaK = max(self.distMat[i, k], self.distMat[k, j])
                    if viaK < abs(self.distMat[i, j]) and abs(self.distMat[i, j])!= 0:
                        self.distMat[i, j] = -viaK
            for i in range(len(self.distMat)):
                for j in range(len(self.distMat)):
                    if self.distMat[i, j] < 0:
                        self.adjMat[i, j] = 0

    def non_gapped(self, root, penalty):

        trans_dict = defaultdict(list)
        final_dict = dict()
        for key in self.edges.keys():
            trans_dict[key[0]].append((key[1], self.edges[key]))
            trans_dict[key[1]].append((key[0], list(map(lambda x: x[::-1], self.edges[key]))))

        visited, queue = set(), deque([root])
        while queue:
            vertex = queue.popleft()
            print('Start: ', vertex)
            save_list = list()
            for neighbour in trans_dict[vertex]:
                if neighbour[0] not in visited:
                    visited.add(neighbour[0])
                    if ''.join(neighbour[1]).count('-') <= penalty and neighbour not in save_list:
                        save_list.append(neighbour)
            if len(save_list) > 0:
                if vertex not in final_dict.keys():
                    final_dict[vertex] = save_list
                else:
                    final_dict[vertex] = list(chain(final_dict[vertex], save_list))
            for neighbour in save_list:
                queue.append(neighbour[0])
        self.edges = final_dict

    def export(self, title):
        export = pd.DataFrame({'Start_vertex': list(map(lambda x: x[0], list(self.edges.items()))), \
                               'End_vertex': list(map(lambda x: x[1][0][0], list(self.edges.items()))), \
                               'Mismatches': list(map(lambda x: ', '.join(x[1][0][1]), list(self.edges.items())))})
        export.to_csv(title, sep='\t', header=True)
        return export

    def visualize(self):
        vis = Graph(comment='Cry toxin amino acid substitution graph', strict=True)
        for start, acc in (self.edges.items()):
            vis.node(start, label=f"{start}")
            for pair in acc:
                end, mismatches = pair[0], pair[1]
                vis.node(end, label=f"{end}")
                vis.edge(start, end, label=", ".join(mismatches))

        vis.view()
        vis.save()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='graph_viualization')
#    parser.add_argument('-ai', help='Name/path of your fasta file', metavar='File',
#                        type=str, required=True)
 #   parser.add_argument('-alg', help='Preferred aligner name and/or absolute path', metavar='Program',
  #                      type=str, required=True)
    parser.add_argument('-ao', help='Name/path to your alignment file', metavar='File',
                        type=str, required=True)
    parser.add_argument('-r', help='Root node', metavar='String', 
                        type=str, required=True)
    parser.add_argument('-o', help='Adjacency list report', metavar='File', 
                        type=str, required=True)
    parser.add_argument('-p', help='Maximum allowed number of indels between adjacent nodes', metavar='Int',
                        type=int, required=True)
    parser.set_defaults(feature=True)

    args = parser.parse_args()
    ao, root, penalty, title = args.ao, args.r, args.p, args.o
 #   subprocess.run(f"{alg} --quiet {ai} > {ao}", shell=True)
    
    g = AAGraph(ao)
    g.upload_vertices()
    g.matBuild()
    g.TR()
    g.non_gapped(root, penalty)
    g.export(title)
    g.visualize()
