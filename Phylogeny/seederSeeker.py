### The output of this script is not a phylogenetic tree. But we use it to explore the phylogenetic relationships between samples
from sys import stdin
import operator
import copy
DEBUG=False

class Variant:
    def __init__(self,id,afs,threshold):
        self.id=id
        self.afs=afs
        self.binary=[1 if i >=threshold else 0 for i in self.afs]
        self.num=sum([self.binary[i]*(2**i) for i in range(len(self.binary))])

class Cluster:
    def __init__(self,id):
        self.id=id
        self.preparent=None
        self.parent=None
        self.children=[]
        self.has_child=False
        self.variants=[]
    def addVariants(self,variant_ids):
        self.variants=variant_ids
    def addMostLikelyParent(self, parent_id,score):
        self.preparent=(parent_id,score)
    def addParent(self,parent_id):
        self.parent=parent_id
    def addChildren(self, children_ids):
        self.children=children_ids

class Node:
    def __init__(self,id):
        self.id=id
        self.pattern=None
        self.map2samples=[]
        self.possible_sample=None

class Sample:
    def __init__(self,id):
        self.id=id
        self.variants={}
        self.best_match_node=None
        self.clusters=[]

def binary2Num(a):
    return sum([a[i]*(2**i) for i in range(len(a))])

def Cluster2Node(v_id,n_id):
    node=Node(n_id)
    v_list=[v_id]
    while(clusters[v_id].preparent!=None):
        v_list.append(clusters[v_id].preparent[0])
        v_id = clusters[v_id].preparent[0]
    node.pattern=[1 if i in v_list else 0 for i in range(num_c) ]
    return node

def mapSample2NodeFunc1():
    for s in samples:
        max_score=-1
        for n in nodes:
            score=sum(s.clusters[i]==n.pattern[i] for i in range(num_c))
            if score> max_score:
                max_score=score
                s.best_match_node=(n.id,max_score)
        nodes[s.best_match_node[0]].map2samples.append(s.id)
def mapSample2Node():
    for s in samples:
        if DEBUG:
            print("~~~~~~~~~~~~~")
            print('sample',s.id)
        max_score=-1
        for n in nodes:
            score=sum([int(s.clusters[i]==n.pattern[i])*len([c for c in clusters if c.id == i][0].variants) for i in range(num_c)])
            if DEBUG:
                print('node', n.id)
                print('pattern', n.pattern)
                print('score', score)
            if score> max_score:
                max_score=score
                s.best_match_node=(n.id,max_score)
        if DEBUG:
            print('best match is node %d and score is %d' % (s.best_match_node[0], s.best_match_node[1]))
        nodes[s.best_match_node[0]].map2samples.append(s.id)
def mapNode2Sample():
    for n in nodes:
        max_score=-1
        for s in samples:
            score=sum(s.clusters[i]==n.pattern[i] for i in range(num_c))
            if score>max_score:
                max_score=score
                n.possible_sample=(s.id,max_score)

# def mapSample2Node():
#     pass


def outPutGraphVis():
    print('digraph{')
    for n in nodes:
        if len(n.map2samples)>0:
            print(str(n.id)+' [label="'+' '.join([samplenames[i] for i in n.map2samples])+'"];')
        else:
            #print(str(n.id)+' [label="'+ 'possible '+samplenames[n.possible_sample[0]]+'"];')
            print(str(n.id)+' [label="inferred ancestor"];')
    for n in nodes:
        if clusters[n.id].preparent!=None:
            print(clusters[n.id].preparent[0],'->',n.id,';')
    print('}')

#parameters:
cutoff=5 # num of variants in each cluster
threshold=0.05
#output='print'
output='graphvis'

# initiation
raw_matrix=[]
row=0
variants=[]
samples=[]
samplenames = []
clusters=[]
variant2num=dict()
v2n_list=[]
kept_v=[]
num_s=0
for line in stdin:
    line=line.strip()
    a=line.split(' ')
    if a[0][0] == "#":
        samplenames = a
    else:
        a=[float(i) for i in a]
        raw_matrix.append(a)
        V=Variant(row,a,threshold)
        variants.append(V)
        v2n_list.append(V.num)
        if V.num not in variant2num:
            variant2num[V.num]=[V.id]
        else:
            variant2num[V.num].append(V.id)
        row+=1

num_s=len(raw_matrix[0])
num_v=len(raw_matrix)
samples=[Sample(i) for i in range(num_s)]
for s in range(num_s):
    for i in range(num_v):
        samples[s].variants[i]=variants[i].afs[s]

# print("samples")
# for s in samples:
#     print("sample id",s.id)
#     print("sample variants", s.variants)

for key in variant2num:
    if len(variant2num[key])>cutoff and key not in kept_v:
        kept_v.append(key)
index = [v2n_list.index(i) for i in kept_v]
matrix=[v.binary for v in variants if v.id in index]
num_c=len(matrix)
for s in range(num_s):
    samples[s].clusters=[matrix[i][s] for i in range(num_c)]

clusters=[Cluster(i) for i in range(len(matrix))]
[c.addVariants(variant2num[binary2Num(matrix[c.id])]) for c in clusters]

for i in range(len(clusters)-1):
    for j in range(i+1, len(clusters)):
        if min([matrix[i][t]-matrix[j][t] for t in range(num_s)]) >=0:
            clusters[i].has_child=True
            score = sum([matrix[i][t]==matrix[j][t] for t in range(num_s)])
            if clusters[j].preparent != None:
                if clusters[j].preparent[1]<score:
                    clusters[j].addMostLikelyParent(i,score)
            else:
                clusters[j].addMostLikelyParent(i,score)
        elif min([matrix[j][t]-matrix[i][t] for t in range(num_s)]) >=0:
            clusters[j].has_child=True
            score = sum([matrix[i][t]==matrix[j][t] for t in range(num_s)])
            if clusters[j].preparent != None:
                if clusters[i].preparent[1]<score:
                    clusters[i].addMostLikelyParent(j,score)
            else:
                clusters[i].addMostLikelyParent(j,score)

# print("clusters")
# for c in clusters:
#     print("cluster id", c.id)
#     print("included variants",c.variants)
#     print("preparent", c.preparent)

nodes=[Cluster2Node(v.id,v.id) for v in clusters]
#print("number of nodes",len(nodes))
# for n in nodes:
#     print("node id",n.id)
#     print("node pattern",n.pattern)
#     print("number of element in pattern",len(n.pattern))
#     print("possible_sample",n.possible_sample)
mapSample2NodeFunc1()
mapNode2Sample()
if output=='graphvis':
    outPutGraphVis()
if output=='print':
    print('~~~~~~~~`cluster tree~~~~~~~')
    for v in clusters:
        print(v.id)
        print(v.preparent)
        print(v.has_child)
    print('number of clusters',num_c)
    for s in samples:
        print(samplenames[s.id])
        print('map to node', s.best_match_node[0])
        print('map score', s.best_match_node[1])
