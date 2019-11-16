import pandas as pd
import numpy as np
from Bioinformatics import Bioinformatics as Bio
import matplotlib.pyplot as plt
import itertools


# Excercise 1_3_10 construct a 4-universal string...
# Bioinformatics Algorthims by Phillip Compeau and Pavel Pavzner, ch. 3

# This example regards making a 'universal string' 
# for example 0001110100 is a 3-universal string because it contains all of
# the 3-base combinations (000, 001, 011, 111, 110, 101, 010, and 100)

bases=['0','1',]
k=3
all_combos=[''.join(p) for p in itertools.product(bases, repeat=k)]


# First we make a relation graph where each room shows the entrance to other rooms with notation
# 101 -> 010,011
# This means 1010 for example is a way where room 101 can lead to room 010
# Here we have an overlap of 2, this whole example is eventually useful for
# Bioinformatics, but I will ignore that for now. 

ans=Bio.get_relation_graph(all_combos)
for relation in ans:
    print(relation)
    
# Get the adjacency matrix that shows the relation between elements that are one charecter apart
# i.e if (000, 001, 011, 111, 110, 101, 010, and 100)=(a,b,c,d,e,f,g,h)
# Adjacency matrix=
#   a  b  c  d  e  f  g  h
# a 0  1  0  0  ect...
# b
# c
# d
# e
# f
# g
# h 1  1 0  0  0  0  0  0
    
    
adj_mat=Bio.get_adjacency_matrix(all_combos)



# Label rooms with letters but use the relation graph
rooms=all_combos
travel_route=[]

history_mat=adj_mat # Matrix that tracks where to go during the alogirthim
# if a 1 is present, check that room, if that room exists on the route
# go back, mark that door with a zero and try the next one
# if out of doors, go back one room on the travel route and reset that row
# to same as adj_mat

# Start in first room, we know this is the start ( or can be so)
# will ignore how to find the start for now

def get_all_doors(room,rooms,adj_mat):
    doors=[]
    # get index of rooms, this will tell you what row of adj_mat to use
    room_index=rooms.index(room)
    adj_row=adj_mat[room_index][:]
    
    for i,possible_room in enumerate(rooms):
        if (adj_row[i]==1):
            doors.append(possible_room)
            
    return doors

#input starting room
travel_route=[rooms[0]]
current_room=rooms[0]
current_room_index=rooms.index(current_room)

history_mat=adj_mat.copy()

try_counter= 0       # Track how many times a route is tried
#print("rooms are")
#print(rooms)
break_condition=1
while(not (len(travel_route)==len(rooms))):
    try_counter=try_counter+1
    previous_travel_route=travel_route
    # Now elimate any doors we cant go through, we need to only travel through each element once
    doors=get_all_doors(current_room,rooms,history_mat)   
    for door in doors:
        # Logic conditions to check
        door_leads_to_room_already_been=(door in travel_route)
        never_been_there_before_condition=not (door in travel_route)
        room_index=rooms.index(door)


        
        if(door_leads_to_room_already_been):
            # get the index for the room so we can
            # erase the choice to go here on the next pass
            history_mat[current_room_index][room_index]=0
            #print('setting '+str(current_room_index)+","+str(room_index)+' to 0')
            #print(history_mat[current_room_index][:])
            #break
            
        #door_check=history_mat[current_room_index][room_index].copy()
        #history_says_ok=(door_check==1)   
            
        if(never_been_there_before_condition):
            current_room=door
            current_room_index=rooms.index(current_room)
            travel_route.append(door)
            #print('ok to go into room '+str(door))
            break
     
    # IF THERE ARE NO MORE POSSIBLE DOORS< GO BACK BRA  
    if ((sum(history_mat[current_room_index][:])==0) ):
        # reset history matrix of current room,
        # remove current room from history matrix
        # of previous room
        previous_room=travel_route[-2]
        previous_room_index=rooms.index(previous_room)
        # remove the previous rooms door to this room
        history_mat[previous_room_index][current_room_index]=0
        # put restore this room...
        history_mat[current_room_index][:]=adj_mat[current_room_index][:].copy()

        back_route=travel_route[0:-1]
        travel_route=back_route
        current_room=previous_room
        current_room_index=rooms.index(current_room)


            ## PRINT CURRENT STATUS

    
travel_route_str=''
for route_room in travel_route:
    travel_route_str=travel_route_str+route_room+','
current_room_options_string=''
for column in history_mat[current_room_index][:]:
    current_room_options_string=current_room_options_string+str(int(column))+','    
print(" route is: ["+travel_route_str+"]")
    
# now combine kmers into a string
ans_string=Bio.get_string_from_ordered_kmers(travel_route)
print('universal string is '+ ans_string)



            
            
            
            
    
            
        
    
   

    




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




    


