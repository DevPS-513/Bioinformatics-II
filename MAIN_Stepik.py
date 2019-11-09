import pandas as pd
import numpy as np
from Bioinformatics import Bioinformatics as Bio
import matplotlib.pyplot as plt
import os
import itertools


# Excercise 1_3_10 construct a 4-universal string...

bases=['0','1',]
k=3
all_combos=[''.join(p) for p in itertools.product(bases, repeat=k)]


in_list=['ATGCG',
'GCATG',
'CATGC',
'AGGCA',
'GGCAT',
'GGCAC']

ans=Bio.get_relation_graph(all_combos)
adj_mat=Bio.get_adjacency_matrix(all_combos)

rooms=all_combos

travel_route=[]

history_mat=adj_mat # Matrix that tracks where to go during the alogirthim
# if a 1 is present, check that room, if that room exists on the route
# go back, mark that door with a zero and try the next one
# if out of doors, go back one room on the travel route and reset that row
# to same as adj_mat

# Start in first room, we know this is the start ( or can be so)
# will ignore how to find the start for now

travel_route=[rooms[0]]


while(len(travel_route)<len(rooms)):
    
    # get index of the current room
    room_index=rooms.index(travel_route[-1])
    print(room_index)
    travel_route.append('000')

    




# Stepik Challenge 1.3 String Spelled by a path
''''
input_file_name='./dataset_198_3.txt'
input_file=open(input_file_name,'r')
mer_list=[]
for i,line in enumerate(input_file.readlines()):
    line=line.rstrip()
    mer_list.append(line)
print(mer_list)
path=Bio.construct_string_from_kmers(mer_list)
print(path)
'''
# Stepik Challenge 1.3.9 print graph relations

'''
input_file_name='./dataset_198_10.txt'
input_file=open(input_file_name,'r')
mer_list=[]


for i,line in enumerate(input_file.readlines()):
    line=line.rstrip()
    mer_list.append(line)


Nmer=len(mer_list)
adj_matrix=np.zeros((Nmer,Nmer))
k=len(mer_list[0])
adj_list=mer_list.copy()
ans_list=[]
ans_list_counter=0
print(mer_list)
for i, merx in enumerate(mer_list):
    prefix_merx=merx[:-1]
    suffix_merx=merx[1:]
    dummy=merx
    adj_list[i] = dummy + ' -> '

    for j,mery in enumerate(mer_list):
        prefix_mery = mery[0:-1]
        suffix_mery = mery[1:]


        #print('merx: '+merx+' sx: '+suffix_merx+' mery: '+mery+'ppy: '+prefix_mery)
        #print(mery)
        # If a match is found add it to the dependency l+merist
        if(prefix_mery==suffix_merx):
            adj_matrix[i][j]=1
            adj_list[i]=adj_list[i]+mery+','

    #remove a comma if it exists

    if (adj_list[i][-1]==','):
        adj_list[i]=adj_list[i][0:-1]

    # If its longer than 'kmer ->'
    # then a relation was added so add it to the answer list
    if (len(adj_list[i])>(k+4)):
        ans_list.append(adj_list[i])


outfile=open('./out_1_3_9_overlap.txt','w')

for item in ans_list:
    outfile.write(item+'\n')


outfile.close
'''




    


