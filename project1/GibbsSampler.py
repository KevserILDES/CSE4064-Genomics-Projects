
import collections
import random

#create profile matrix
def profile_matrix(motif_matrix, removen_pos):
    motif_matrix.remove(motif_matrix[removen_pos])
    # Create profile matrix
    profile = collections.defaultdict(list)
    column = []
    for i in range(0, len(motif_matrix[0])):
        col = ""
        for k in range(0, len(motif_matrix)):
            col += motif_matrix[k][i]
        column.append(col)
    for i in range(0, len(motif_matrix[0])):
        for letter in ['A', 'C', 'G', 'T']:
            letter_count = column[i].count(letter) + 1  # add 1 to avoid zero results
            profile[letter].append(letter_count / (len(motif_matrix[0]) * 2))

    return profile

#calculate probability from profile matrix
def calculate_prob(dna, removen_pos, k_value, profile):
    prob_list = []
    for index in range(0, len(dna[removen_pos]) - k_value + 1):
        prob = 1
        i = 0
        for letter in dna[removen_pos][index:index + k_value]:
            prob *= profile[letter][i]
            i += 1
        if prob > 0:
            prob_list.append(prob)

    return prob_list


#get one motif from probability list
def choosen_motif_from_prob(probs, dna, k):
    sums = float(sum(probs))
    random_prob = random.random()
    my_sum = 0.0
    for i in range(len(probs)):
        my_sum += probs[i]
        if random_prob <= my_sum / sums:
            return dna[i:i + k]


#calculate score
def score(motif_matrix):
    s = []
    column = []
    # get columns from motif matrix
    for i in range(0, len(motif_matrix[0])):
        col = ""
        for k in range(0, len(motif_matrix)):
            col += motif_matrix[k][i]
        column.append(col)
    # find most freq letter from columns and calculate score
    for i in column:
        max_letter = 0
        for letter in ['A', 'C', 'G', 'T']:
            letter_count = i.count(letter)
            if letter_count > max_letter:
                max_letter = letter_count
        s.append(len(i) - max_letter)
    return sum(s)


# run gibbs algorithm
def gibbs_sampler(dna, k, iter_num):
    position = [random.randint(0, len(x) - k) for x in dna] #get random position from lines for create motirf matrix
    i = 0
    motif_matrix = []
    #create motif matrix
    for line in dna:
        motif_matrix.append(line[position[i]:position[i] + k])
        i += 1
    # at the beginning assign initial motifs to best motifs
    bests = motif_matrix.copy()
    # get initial score
    last_score = score(bests)
    while iter_num < 50:
        #remove randomly one motif from motirf matrix
        removen_pos = random.randint(0, len(motif_matrix) - 1)
        profile = profile_matrix(motif_matrix, removen_pos)
        probs = calculate_prob(dna, removen_pos, k, profile)
        #replace the deleted motif with the motif with best score
        choosen_motif = choosen_motif_from_prob(probs, dna[removen_pos], k)
        motif_matrix.insert(removen_pos, choosen_motif)
        new_score = score(motif_matrix)
        if new_score < last_score:
            bests = motif_matrix.copy()
            last_score = new_score
            iter_num = 0
        else:
            iter_num += 1

    return bests



#run gibbs algorithm 50 times and get best result
def run_gibbs(dna, k, iter_num):
    motifs = gibbs_sampler(dna, k, iter_num)
    sum_scores = 0
    max_score = 0
    for i in range(50):
        newmotifs = gibbs_sampler(dna, k, iter_num)
        sum_scores += score(newmotifs)
        if score(newmotifs) > max_score:
            max_score = score(newmotifs)
        if score(newmotifs) < score(motifs):
            motifs = newmotifs.copy()
    return motifs, score(motifs), max_score, sum_scores/50


#get user input
def get_input():
    while True:
        try:
            k_value = int(input("Please enter value of k (9,10 or 11): "))
        except ValueError:
            print("Wrong input")
            continue
        if not (k_value == 9 or k_value == 10 or k_value == 11):
            print("Values of k can be 9, 10 or 11.")
        else:
            file_name = input("Please enter txt file name: ")
            file = open(file_name, 'r')
            dna = [line.strip() for line in file]
            file.close()
            # Make sure the sequences are all upper case
            dna = [seq.upper() for seq in dna]
            return dna, k_value


dna, k_value = get_input()
iter_num = 0

#print results
motifs, score, max_score, avg_score = run_gibbs(dna, k_value, iter_num)
print("Motifs", "\n--------------")
for i in motifs:
    print(i)
print("\nScore Result","\n--------------")
print("Best Score : ", score, "\nMax Score: ", max_score, "\nAverage Score: ", avg_score)



