import pandas as pd
import numpy as np
from Bioinformatics import Bioinformatics as Bio
import matplotlib.pyplot as plt

#This Code calculates the Median String Example from
# http://www.csbio.unc.edu/mcmillan/Comp555S16/Lecture05.html
# Where the Motif TAGATCCGAA exists in 10 strands of DNA
# The motif is noisy, meaning it can appear slightly different in each strand
# i.e in strand 2 it appears as TGGATCCGAA but is considered a motif of TAGATCCGAA

print('The Example at http://www.csbio.unc.edu/mcmillan/Comp555S16/Lecture05.html ' )
print('is looking for embedded noisy motif TAGATCCGAA from:')
print()

infile=open('.\DNA_strands.txt',"r")
DNA_strands_o=infile.readlines()
DNA_strands=['']*len(DNA_strands_o)
for i,strand in enumerate(DNA_strands_o):
    new_strand=strand.rstrip()          # Remove the new line charecter
    DNA_strands[i]=new_strand.upper()   # Convert to Uppercase if needed
    print(DNA_strands[i])

print()
print('This code will compare results from three different methods to find TAGATCCGAA ')
# Question asks for RandomizedMotifSearch
k=10
# Method 1, Median String
print('Method\t\tConsensus\tEntropy')
print('------------------------------------------------------')
winning_list=Bio.get_MedianString_motifs(k,DNA_strands,0)
profile_matrix=Bio.get_ACGT_profile_matrix_from_motifs_psudocount(k,winning_list)
entropy=Bio.get_entropy_from_profile_matrix_ACTG(k,profile_matrix)
consensus_mer=Bio.get_consensus_mer(winning_list,profile_matrix)
print('MediaString\t'+str(consensus_mer)+'\t'+str(round(entropy,4)))



iterations=100
limit_iterations=100  # limit number of times it sub-iterates in case it gets stuck...
motifs_and_score=Bio.Randomized_Motif_Search(iterations,limit_iterations,k,DNA_strands)
profile_matrix=Bio.get_ACGT_profile_matrix_from_motifs_psudocount(k,motifs_and_score[0])
entropy=Bio.get_entropy_from_profile_matrix_ACTG(k,profile_matrix)
consensus_mer=Bio.get_consensus_mer(winning_list,profile_matrix)
print('Randomized\t'+str(consensus_mer)+'\t'+str(round(entropy,4)))


iterations=100
limit_iterations=100  # limit number of times it sub-iterates in case it gets stuck...
motifs_and_score=Bio.Gibbs_Sampler(k,DNA_strands,iterations,limit_iterations)
profile_matrix=Bio.get_ACGT_profile_matrix_from_motifs_psudocount(k,motifs_and_score[0])
winning_list=motifs_and_score[0]
consensus_mer=Bio.get_consensus_mer(winning_list,profile_matrix)
print('Gibbs Sampler\t'+str(consensus_mer)+'\t'+str(round(motifs_and_score[1],4)))






    

