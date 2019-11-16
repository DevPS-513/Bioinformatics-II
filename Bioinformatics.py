
# A list of functions used throught the Bionformatics course
# on courser.org
class Bioinformatics:
    
       
    # This allows the use of Bioinformatics().HammingDistance(...)
    # So that within a method you can call another method
    # This was done so that a helpful drop down list appears when you type
    # 'Bioinformactics().' in the IDE i'm using..


    def __new__(self): 
        return self
    

    
    def HammingDistance(text1,text2):
        # Get the Hamming Distance between two equally sized strings
        # i.e 'ACG' and 'AAA', d=2
        # 'AAA', and 'AAA', d=0
        hamming_counter=0
        for j,char in enumerate(text1):
            if (char!=text2[j]):
                hamming_counter=hamming_counter+1    
    
        return hamming_counter
    
  
    def get_all_kmers(k:int, DNA_list:list)->list:
        # This function returns all k-mers within a list of DNA strings
        # For example (k=3,DNA_list=['AAATTGACGCAT'])
        # will return ['AAA', 'AAT', 'ATT', 'TTG', 'TGA', 'GAC', 'ACG', 'CGC', 'GCA', 'CAT']        
        # first if only one string is passed, convert to a list
        if(type(DNA_list)==str):
            DNA_list=[DNA_list]
            
        N_bases=len(DNA_list[0])
        list_of_mers=[]
        for strand in DNA_list:
            N_mers=N_bases-k+1 # kmers, imagine scanning from start to the last one, only the start of the last scan is included
            for j in range(0,N_mers):    
                list_of_mers.append(strand[j:j+k])    
        return list_of_mers

    def get_string_from_ordered_kmers(mer_list):
        # Assumeing a list is given in sequence
        # This will stick together each string with only one charecter
        # that does not overlap. For example [AAA+AAT=AAAT]
        string=''
        for mer in mer_list:
            string=string+mer[0]

        string=string+mer_list[-1][1:]

        return string
   
    def get_all_composition_kmers(k: int, DNA_list: list) -> list:
        # This function returns all k-mers within a list of DNA strings
        # For example (k=3,DNA_list=['AAATTGACGCAT'])
        # will return ['AAA', 'AAT', 'ATT', 'TTG', 'TGA', 'GAC', 'ACG', 'CGC', 'GCA', 'CAT']

        # Sometimes only one DNA string is passed in this case convert to a list...
        if (type(DNA_list) == str):
            DNA_list = [DNA_list]

        N_bases = len(DNA_list[0])
        # print("length of DNA is "+str(N_DNA))
        list_of_mers = []
        for strand in DNA_list:
            N_mers = N_bases/k  # kmers, imagine scanning from start to the last one, only the start of the last scan is included
            N_mers_int=round(N_bases/k)
            if((N_mers-N_mers_int)>0):
                print('Does not divide evenly into k-mers!')

            N_mers=int(N_mers)
            for j in range(0, N_mers):
                print(j)
                list_of_mers.append(strand[(j-1)*k:(j-1)*k+k])
        return list_of_mers
    
    def get_relation_graph(mer_list):
        # This code displays any relations between a list of kmers,i.e 
        # ATGCG,GCATG,CATGC,AGGCA,GGCAT,GGCAC
        # would output as strings.
        #CATGC -> ATGCG
        #GCATG -> CATGC
        #GGCAT -> GCATG
        #AGGCA -> GGCAC,GGCAT
        import numpy as np
        

        Nmer=len(mer_list)
        adj_matrix=np.zeros((Nmer,Nmer))
        k=len(mer_list[0])
        adj_list=mer_list.copy()
        ans_list=[]
        ans_list_counter=0
        
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
                c1=(prefix_mery==suffix_merx)
                c2=(merx!=mery) # dont include self refer
                if(c1 and c2):
                    adj_matrix[i][j]=1
                    adj_list[i]=adj_list[i]+mery+','
        
            #remove a comma if it exists
        
            if (adj_list[i][-1]==','):
                adj_list[i]=adj_list[i][0:-1]
        
            # If its longer than 'kmer ->'
            # then a relation was added so add it to the answer list
            if (len(adj_list[i])>(k+4)):
                ans_list.append(adj_list[i])

        return ans_list
        
        
        
    def get_adjacency_matrix(mer_list):
        # This code displays any relations between a list of kmers,i.e 
        # ATGCG,GCATG,CATGC,AGGCA,GGCAT,GGCAC
        # would output as strings.
        #CATGC -> ATGCG
        #GCATG -> CATGC
        #GGCAT -> GCATG
        #AGGCA -> GGCAC,GGCAT
        import numpy as np


        Nmer=len(mer_list)
        adj_matrix=np.zeros((Nmer,Nmer))
        k=len(mer_list[0])
        adj_list=mer_list.copy()
        ans_list=[]
        ans_list_counter=0

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
                c1=(prefix_mery==suffix_merx)
                c2=(merx!=mery) # dont include self refer
                if(c1 and c2):
                    adj_matrix[i][j]=1
                    adj_list[i]=adj_list[i]+mery+','

            #remove a comma if it exists

            if (adj_list[i][-1]==','):
                adj_list[i]=adj_list[i][0:-1]

        # If its longer than 'kmer ->'
        # then a relation was added so add it to the answer list
        if (len(adj_list[i])>(k+4)):
            ans_list.append(adj_list[i])

        return adj_matrix
    
    
    
    def get_profile_most_probable_kmer_from_list(k:int,strand:str,profile_matrix:list)->str:    
        # This function uses a profile matrix, which represents the probability
        # of each DNA base 'A' 'C' 'T' 'G' within a kmer of k length
        # it extracts every kmer in the DNA strand and calculates its probability
        # Format of probability matrix is assumed to be ACGT, for each position in the kmer
        # i.e if probability_matrix=[[0.25,0.25,0.5,0],[0.125,0.125,0.5,0.25]]
        # for kmer='AT', k=2, prob=profile_matrix[1][1]*profile_matrix[2][3]
        mer_scores=[]
        list_of_mers=Bioinformatics().get_all_kmers(k,[strand])
        prob_mat=[0]*k
        prob_occur=1
        
        for mer in list_of_mers:
    
            prob_occur = 1
            for j, char in enumerate(mer):
                if (char == "A"):
                    prob_mat[j]=profile_matrix[j][0]
                if (char == "C"):
                    prob_mat[j]=profile_matrix[j][1]
                if (char == "G"):
                    prob_mat[j]=profile_matrix[j][2]
                if (char == "T"):
                    prob_mat[j]=profile_matrix[j][3]
    
            for p in prob_mat:
                prob_occur=prob_occur*p
    
            #Track the probability score for that mer
            mer_scores.append(prob_occur)
    
        # find maximum mer
        max_mer_score=max(mer_scores)
        max_mer_score_loc=mer_scores.index(max_mer_score)
        max_mer=list_of_mers[max_mer_score_loc]
        
        return max_mer
        
    
    def get_ACGT_profile_matrix_from_motifs_psudocount(k:int,list_of_motifs:list)->list:
       
        # This returns a profile matrix from a list of motifs, i.e
        # list_of_motifs =['AA','AC']
        # count_matrix=[[2,0,0,0],[1,1,0,0]]
        # Profile_matrix=[[1,0,0,0],[0.5,0.5,0,0]]
        # Psudocounts will pretend like each letter has at least one occureance
        # count_matrix=[[2,0,0,0],[1,1,0,0]]
        # count_matrix_psudo=[[3,1,1,1],[2,2,1,1]]
        Nmers=len(list_of_motifs)
    
        # This simply initalizes the profuile and count matrices
        # likely a simpler way to do this but this will work
        profile_matrix=[0]*k   
        count_matrix=[0]*k  
        for i in range(k):
            profile_matrix[i]=[0,0,0,0] # Add four ACGT values to each profile column
            count_matrix[i]=[0,0,0,0] # Add four ACGT values to each profile column
        # Start finding the profile_matrix        
        for i,mer in enumerate(list_of_motifs):
            for j,char in enumerate(mer):
                if(char=="A"):
                    count_matrix[j][0]=count_matrix[j][0]+1
                if(char=="C"):
                    count_matrix[j][1]=count_matrix[j][1]+1
    
                if(char=="G"):
                    count_matrix[j][2]=count_matrix[j][2]+1
    
                if(char=="T"):
                    count_matrix[j][3]=count_matrix[j][3]+1 
         
    # Now convert to probability
    # Because we are adding one to each position, its as if there was twice as many mers
    # with the second list being evenly distributed
        for q in range(k):
            profile_matrix[q][0]=(count_matrix[q][0]+1)/(Nmers*2)
            profile_matrix[q][1]=(count_matrix[q][1]+1)/(Nmers*2)
            profile_matrix[q][2]=(count_matrix[q][2]+1)/(Nmers*2)
            profile_matrix[q][3]=(count_matrix[q][3]+1)/(Nmers*2)       
    
    
        return profile_matrix      
    
    def get_ACGT_profile_matrix_from_motifs(k,list_of_motifs):
        # This returns a profile matrix from a list of motifs, i.e
        # list_of_motifs =['AA','AC']
        # count_matrix=[[2,0,0,0],[1,1,0,0]]
        # Profile_matrix=[[1,0,0,0],[0.5,0.5,0,0]]

       
        Nmers=len(list_of_motifs)
        profile_matrix=[0]*k   
        count_matrix=[0]*k  
        for i in range(k):
            profile_matrix[i]=[0,0,0,0] # Add four ACGT values to each profile column
            count_matrix[i]=[0,0,0,0] # Add four ACGT values to each profile column
           
        for m,mer in enumerate(list_of_motifs):
            for j,char in enumerate(mer):

                if(char=="A"):
                    count_matrix[j][0]=count_matrix[j][0]+1
                if(char=="C"):
                    count_matrix[j][1]=count_matrix[j][1]+1
                if(char=="G"):
                    count_matrix[j][2]=count_matrix[j][2]+1
                if(char=="T"):
                    count_matrix[j][3]=count_matrix[j][3]+1
         
        # Now convert to probability
        for q in range(k):
            profile_matrix[q][0]=(count_matrix[q][0])/Nmers
            profile_matrix[q][1]=(count_matrix[q][1])/Nmers
            profile_matrix[q][2]=(count_matrix[q][2])/Nmers
            profile_matrix[q][3]=(count_matrix[q][3])/Nmers       
    
    
        return profile_matrix  

    def get_entropy_from_ACTG(probs:list):
        entropy=0
        import math
        for p in probs:
            if(p==0):
                entropy=entropy+0
            else:
                entropy=entropy+p*math.log2(p)
                
        return -entropy


    def get_entropy_from_profile_matrix_ACTG(k,profile_matrix):
        # this will calcualate the total entropy of a profile matrix
        #i.e profile_matrix=[[0.25,0.25,0.25,0.25],[0,0.5,0.5,0]]
        # e1=entropy(profile_matrix[0]=sum(profile_matrix[0][:]*log_2(profile_matrix[0][:]))
        # e2=entropy(profile_matrix[1]=sum(profile_matrix[1][:]*log_2(profile_matrix[1][:]))
        # entropy=e1+e2   
   
        entropy=0 
        for col in range(k):
            entropy=entropy+Bioinformatics().get_entropy_from_ACTG(profile_matrix[col][:])  
        return entropy

    def get_list_of_kmer_probabilities(k,profile_matrix,DNA_strand):
        # Now given a DNA strand, seperate into all possible kmers
        # and get the probability each kmer occurs given a profile_matrix
        list_of_kmers=Bioinformatics().get_all_kmers(k,DNA_strand)
        Nmers=len(list_of_kmers)
    
        list_of_probabilities=[1]*Nmers
        sum_probs=0
        for j,mer in enumerate(list_of_kmers):
            
            for i,char in enumerate(mer):
                if (char=='A'):
                    list_of_probabilities[j]=list_of_probabilities[j]*profile_matrix[i][0]
                if (char=='C'):
                    list_of_probabilities[j]=list_of_probabilities[j]*profile_matrix[i][1]
                if (char=='G'):
                    list_of_probabilities[j]=list_of_probabilities[j]*profile_matrix[i][2]
                if (char=='T'):
                    list_of_probabilities[j]=list_of_probabilities[j]*profile_matrix[i][3]
                    
            sum_probs=sum_probs+list_of_probabilities[j]
        
        for q,prob in enumerate(list_of_probabilities):
            list_of_probabilities[q]=list_of_probabilities[q]/sum_probs
        
        return list_of_probabilities


    def get_consensus_mer(list_of_kmers,profile_matrix): 
        # Given a profile matrix return the most likely kmer
        # for example profile_matrix=[[0.6,0.4,0,0],[0,0,1,0]]
        # consensus_mer='AG'
        # Just the first one is given if there is more than one possibility.
        k=len(list_of_kmers[0])
        
        concensus_string=''   
        # Now go through each column and generate letter where maximum is...
        base_list=['A','C','G','T']
        for i in range(k):
            prob_list=profile_matrix[i][:]
            max_loc=prob_list.index(max(prob_list))
            max_base=base_list[max_loc]
            
            concensus_string=concensus_string+max_base
                
        return concensus_string
    



    def generate_all_kmers(k:int):
        # Generate all possible kmers, i.e 4^k possibilities
        # k=2, list=['AA','AC','AG','AT','CA','CC','CG','CT',...ect]
        import itertools
        bases=['A','C','G','T']
        all_kmers=[''.join(p) for p in itertools.product(bases, repeat=k)]
        
        return all_kmers

    def get_minimum_hamming_mer_and_score(kmer: str,DNA:str)->tuple:
        # looks at DNA string and finds the minimum hamming distance
        # of kmer with all possible kmers of the same length in DNA
    
        k=len(kmer)
        N_bases=len(DNA)
      
        N_mers=N_bases-k+1 # kmers, imagine scanning from start to the last one, only the start of the last scan is included
        list_of_mers=[]
        list_of_hamming_distances=[]
        for j in range(N_mers):    
            list_of_mers.append(DNA[j:j+k])
            current_mer=list_of_mers[j]
            list_of_hamming_distances.append(Bioinformatics().HammingDistance(kmer,current_mer))
            
    
        min_mer_loc=list_of_hamming_distances.index(min(list_of_hamming_distances))        
        
        min_mer=list_of_mers[min_mer_loc]
        min_hamming=list_of_hamming_distances[min_mer_loc]
           
        return (min_mer,min_hamming)
     
    
    def get_MedianString_motifs(k:int, DNA_strands:list, brute:int)->tuple:
        # Will first assume all patterns to try are generated from
        # THe DNA strings only...
        
            
        # Now we will go through every DNA strand and find the minimum
        # Hamming Distance in each strand, the sum will be the score
        # Hamming Distance is the number of chars that ar different
        # i.E if kmer is 'AAA' 
        #   TTACCTT[AAC]   =1
        #   G[ATA]TCTGTG   =1
        #   [ACG]GCGTTCG   =2
        #   CCCTAA[AGA]G   =1
        # score=5 for 'AAA'
        if (brute==0):
            list_of_all_kmers=Bioinformatics().get_all_kmers(k,DNA_strands)
        else:
            list_of_all_kmers=Bioinformatics().generate_all_kmers(k)
    
    
        Nstrands=len(DNA_strands)
        Nmers=len(list_of_all_kmers)
        best_score=Nstrands*Nmers # Theoretically should be impossible to be worse then this score
        
        for q,kmer in enumerate(list_of_all_kmers):
            score=0
            motif_list=[]
            for strand in DNA_strands:
                mer_and_score=Bioinformatics().get_minimum_hamming_mer_and_score(kmer,strand)
                score=score+mer_and_score[1]
                motif_list.append(mer_and_score[0])
            
            # After summing over all strands check the score of mer
            if(score<best_score):
                best_score=score
                profile_matrix=Bioinformatics().get_ACGT_profile_matrix_from_motifs_psudocount(k,motif_list)
                ConsensusMer=Bioinformatics().get_consensus_mer(motif_list,profile_matrix)
                winning_list=motif_list
                
                
                
        return winning_list
    
    def Randomized_Motif_Search(iterations:int,limit_iterations:int,k:int,DNA_strands:list):
        import random 

        Nt=len(DNA_strands)
        score_track=[]        
        iterate_flag=1  
        num_of_kmers_per_strand=len(DNA_strands[0])-k+1
        
        list_of_all_kmers=['']*Nt  
        for j,strand in enumerate(DNA_strands):
            kmers_in_strand=Bioinformatics().get_all_kmers(k,strand)
            list_of_all_kmers[j]=kmers_in_strand # Add four ACGT values to each profile column
        

        for q in range(iterations):
        # for each iterartion create a new random start.
            current_motifs=['']*Nt 
            for i in range(Nt):
                current_motifs[i]=list_of_all_kmers[i][random.randint(0,num_of_kmers_per_strand-1)]
                    
                    
                
            profile_matrix=Bioinformatics().get_ACGT_profile_matrix_from_motifs_psudocount(k,current_motifs)
        
            if(q==0):
                profile_matrix=Bioinformatics().get_ACGT_profile_matrix_from_motifs_psudocount(k,current_motifs)
                entropy=Bioinformatics().get_entropy_from_profile_matrix_ACTG(k,profile_matrix)
                score=entropy
                best_score=score
                
            
            sub_iterations=0
            iterate_flag=1
            
            while((iterate_flag==1) and (sub_iterations<limit_iterations)):
                
                sub_iterations=sub_iterations+1
                current_motifs=['']*Nt
                
                for j,strand in enumerate(DNA_strands):
                    new_mer=Bioinformatics().get_profile_most_probable_kmer_from_list(k,strand,profile_matrix)
                    current_motifs[j]=new_mer
                    
                # Redfine profile matrix after all new mers are added
                profile_matrix=Bioinformatics().get_ACGT_profile_matrix_from_motifs_psudocount(k,current_motifs)
                
                consensus_mer=Bioinformatics.get_consensus_mer(current_motifs,profile_matrix)
                #score=score_motif_list_with_DNA(current_motifs,DNA_strands)
                
                entropy=Bioinformatics().get_entropy_from_profile_matrix_ACTG(k,profile_matrix)
                score=entropy
                #print('score is now '+str(round(score,5))+ ' prev motif[0] '+prev_motifs[0]+' new motif[0] '+current_motifs[0] )
                
                if(score<best_score):
                    winning_motifs=current_motifs
                    winning_mer=consensus_mer
                    best_score=score
                else:
                    iterate_flag=0
                    score_track.append(entropy)

        return (winning_motifs,best_score)
    
    
    def Gibbs_Sampler(k:int,DNA_strands:list, iterations:int,limit_iterations:int):
        import random
        import math
        import numpy as np
        
        Nt=len(DNA_strands)
        num_of_kmers_per_strand=len(DNA_strands[0])-k+1
        
        for q in range(iterations):
        # for each iterartion create a new random start.
            current_motifs=['']*Nt 
            
            # STARTING MOTIFS< PROFILE< AND SCORE
            for i,strand in enumerate(DNA_strands):
                list_of_kmers_from_strand=Bioinformatics.get_all_kmers(k,strand)
                rand_mer_index=random.randint(0,num_of_kmers_per_strand-1)
                current_motifs[i]=list_of_kmers_from_strand[rand_mer_index]
                
          
            if(q==0):
                # Start with some outrageous score that will always be beat                
                best_score=math.pow(k*Nt*num_of_kmers_per_strand,2)
                # STARTING MOTIFS< PROFILE< AND SCORE
           
        
            # Reset sub_iterations counter and the flag    
            sub_iterations=0    # Counts to limit iterations in case there is a problem
            iterate_flag=1      # set to 0 when convergence is found
            
            while((iterate_flag==1) and (sub_iterations<limit_iterations)):
                
                sub_iterations=sub_iterations+1
                
                # Now we replace just ONE mer at a time
                random_strand=random.randint(0,Nt-1)
                reduced_motifs=[] # reset reduced motifs
                # Go through each motif, if not part of the deleted strand then
                # add it to reduced motifs list
                for j in range(0,Nt):
                    if (j!=random_strand):
                        reduced_motifs.append(current_motifs[j])
      
                profile_matrix=Bioinformatics.get_ACGT_profile_matrix_from_motifs_psudocount(k,reduced_motifs)
                list_of_kmer_probs=Bioinformatics.get_list_of_kmer_probabilities(k,profile_matrix,[DNA_strands[random_strand]])
                list_of_hidden_kmers=Bioinformatics.get_all_kmers(k,DNA_strands[random_strand])
         
                Nprobs=len(list_of_kmer_probs)
                new_mer_index=np.random.choice(Nprobs,1,p=list_of_kmer_probs)
                new_mer_index=new_mer_index[0]
                new_mer=list_of_hidden_kmers[new_mer_index]
                
                # SCORE THE NEW MOTIFS THAT JUST HAS ONE NEW MER
                current_motifs[random_strand]=new_mer              
                profile_matrix=Bioinformatics.get_ACGT_profile_matrix_from_motifs_psudocount(k,current_motifs)
        
                #score=score_motif_list_with_DNA(current_motifs,DNA_strands)
                entropy=Bioinformatics.get_entropy_from_profile_matrix_ACTG(k,profile_matrix)
                score=entropy
    
                # Should brake out of loop if guesses do not get better
                if(score<best_score):
                    winning_motifs=current_motifs
                    best_score=score
                    winning_profile=profile_matrix
                else:
                    iterate_flag=1
                    
                
                
        return (winning_motifs,best_score)
