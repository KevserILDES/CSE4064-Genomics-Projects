
import random

"""randomized motif search method takes 3 parameters
dna: our input motifs formatted as lines and all uppercase
k: length of consenuss string (9, 10 or 11)
t: how many iterations ? (given as 50 for this project)
"""
def RandomizedMotifSearch(dna,k,t):
    #create initial motifs for search
    motifs=create_initial_motifs(dna,k)
    #assign initial best motifs and calculate best score
    bestMotifs=motifs.copy()
    bestScore=calculate_score(motifs)
    iterNum=0
    #loop till iterations number reached
    while True:
        #create profile matrix and using this profile create new motifs and calculate new score
        profile=profile_matrix(motifs)
        motifs=create_motifs(profile,dna,k)
        score=calculate_score(motifs)
        #if new score is better update best motifs and score
        if score<bestScore:
            iterNum=0
            bestMotifs=motifs
            bestScore=score
        #increment counter
        else:
            iterNum+=1
        if iterNum==t:
            break
    #return best motifs and score found
    return bestMotifs,bestScore
    
""" Method for creating profile matrix, only one argument which is motifs"""
def profile_matrix(motifs):
    #initialize matrix
    matrix=[[0]*4]*len(motifs[0])
    #genes all with some value, to count occurences
    genes={'A':0,'C':0,'G':0,'T':0}
    #loop through motifs and count occurence of each, first loop horizontal second is vertical
    for i in range(0,len(motifs[0])):
        for j in range(0,len(motifs)):
            #increment corresponding gene's value
            genes[motifs[j][i]]+=1
        #calculate and assign matrix (columns as rows for now)
        matrix[i]=[value/len(motifs) for value in genes.values()]
        #reset values
        genes={'A':0,'C':0,'G':0,'T':0}
    #return profile matrix (which is traverse of created matrix (rows - columns))
    return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

"""Method for calculating probabilities of motifs, 2 arguments which are profile matrix created and motif"""
def calculate_prob(profile,motif):
    prob=1
    #indices of genes in profile matrix
    indices={'A':0,'C':1,'G':2,'T':3}
    for i in range(len(motif)):
        #multiply probability with corresponding value in profile matrix
        prob*=profile[indices[motif[i]]][i]
    #return probability of motif
    return prob

"""Method for calculating score of motifs"""
def calculate_score(motifs):
    #initializations
    score=0
    genes={'A':0,'C':0,'G':0,'T':0}
    for i in range(0,len(motifs[0])):
        for j in range(0,len(motifs)):
            genes[motifs[j][i]]+=1
        #increment score by (10-max occurence of genes)
        score+=(len(motifs)-max(genes.values()))
        genes={'A':0,'C':0,'G':0,'T':0}
    #return score
    return score
    
"""Method for creating motifs,3 arguments:
profile matrix
dna and k value (length of consensus string)
"""
def create_motifs(profile,Dna,k):
    #assignments
    new_motifs=[]
    max_prob=0
    new_motif_of_line=''
    #for each line in dna
    for i in Dna:
        #choose the motif in each line that has max probability
        for j in range(0, len(i) - k):
            #choose motif in length of k till reach end of line by shifting one by one 
            motif=i[j:j+k]
            #calculate probability of motif
            prob=calculate_prob(profile,motif)
            #if new probability is better update new motif
            if(max_prob<prob):
                max_prob=prob
                new_motif_of_line=motif
        #append it to new motifs and reset values for next line
        new_motifs.append(new_motif_of_line)
        max_prob=0
        new_motif_of_line=''
    #return new motifs
    return new_motifs

"""Method for creating initial motifs,2 arguments: dna and length of consensus string"""
def create_initial_motifs(dna,k):
    #initialize
    initial_motifs=[]
    #for each line take a random number and select motif starting from rand with length k
    for i in dna:
        rand=random.randint(0, len(i) - k)
        initial_motifs.append(i[rand:rand+k])
    #return created motifs
    return initial_motifs  
    
#take inputs from user
k_value = int(input("Please enter value of k (9,10 or 11): "))
file_name = input("Please enter txt file name: ")
file = open(file_name, 'r')
#format dna as lines and make all uppercase
dna = [line.strip() for line in file]
dna = [seq.upper() for seq in dna]

#run randomized motif search for dna
bestMotifs=[]
bestScore=9999
maxScore=0
sumScores=0
#run randomized motif search 50 times and get best result
for i in range(50):
    motifs, score = RandomizedMotifSearch(dna,k_value,5)
    sumScores+=score
    if score>maxScore:
        maxScore=score
    if score<bestScore:
        bestMotifs=motifs
        bestScore=score
averageScore=sumScores/50
#print results
for i in bestMotifs:
    print(i)
print(bestScore)
print("Max score: ",maxScore,"\nAverage: ",averageScore)
