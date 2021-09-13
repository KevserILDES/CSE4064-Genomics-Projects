# -*- coding: utf-8 -*-
"""
Kevser İLDEŞ 150116048
Melisa DÖNMEZ - 150116030

"""
import os

global INF

INF = -float("inf")

# default gap values
gap_opening = -1
gap_extension = -0.5
match = 2
mismatch = -1


class SequenceAlignment:
    # assign class values
    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.mismatch_score = mismatch
        self.match_score = match
        self.gap_opn_penalty = gap_opening
        self.gap_ext_penalty = gap_extension

        self.middle_matrix = self.createMatrix(len(seq1) + 1, len(seq2) + 1)
        self.x_matrix = self.createMatrix(len(seq1) + 1, len(seq2) + 1)
        self.y_matrix = self.createMatrix(len(seq1) + 1, len(seq2) + 1)

        self.backtrack_indices = [[[0 for i in range(2)] for j in range(len(seq2) + 1)] for k in range(len(seq1) + 1)]

    # create 0 matrix with given row and column number
    def createMatrix(self, x, y):
        return [[0.0 for i in range(y)] for j in range(x)]

    # assign value to matrix elements
    def initMatrices(self):
        for i in range(1, len(self.seq1) + 1):
            self.middle_matrix[i][0] = INF
            self.x_matrix[i][0] = INF
            self.y_matrix[i][0] = self.gap_opn_penalty + i * self.gap_ext_penalty
        for i in range(1, len(self.seq2) + 1):
            self.middle_matrix[0][i] = INF
            self.x_matrix[0][i] = self.gap_opn_penalty + i * self.gap_ext_penalty
            self.y_matrix[0][i] = INF

    #main alignment method to construct matrices
    def runAlignment(self):
        #initialize all 3 matrices
        self.initMatrices()

        #loop through sequences
        for i in range(1, len(self.seq1) + 1):
            for j in range(1, len(self.seq2) + 1):

                #for upper matrix assignment
                self.x_matrix[i][j] = max(self.x_matrix[i][j - 1] + self.gap_ext_penalty,
                                          self.middle_matrix[i][j - 1] + self.gap_opn_penalty)

                #for lower matrix assignment
                self.y_matrix[i][j] = max(self.y_matrix[i - 1][j] + self.gap_ext_penalty,
                                          self.middle_matrix[i - 1][j] + self.gap_opn_penalty)

                #for middle matrix (main matrix) assignment; first find the matching score and calculate diagonal movement then find max of all alignments
                middle_score = self.calcMatchScore(i - 1, j - 1) + self.middle_matrix[i - 1][j - 1]
                self.middle_matrix[i][j] = max(middle_score,
                                               self.x_matrix[i][j],
                                               self.y_matrix[i][j])

                #according to result store the node from which it come for later backtracking
                if self.middle_matrix[i][j] == self.y_matrix[i][j]:
                    self.backtrack_indices[i][j][0] = i - 1
                    self.backtrack_indices[i][j][1] = j
                elif self.middle_matrix[i][j] == middle_score:
                    self.backtrack_indices[i][j][0] = i - 1
                    self.backtrack_indices[i][j][1] = j - 1
                else:
                    self.backtrack_indices[i][j][0] = i
                    self.backtrack_indices[i][j][1] = j - 1

        #finally, after all matrices filled and backtracking nodes obtained construct alignment and print results
        self.constructAlignment()

    # Calculating the score of matching
    def calcMatchScore(self, row_index, col_index):
        score = 0
        if row_index == 0 or col_index == 0:
            score = 0
        elif self.seq1[row_index] == self.seq2[col_index]:
            score = self.match_score
        else:
            score = self.mismatch_score
        return score

    # traceback
    def constructAlignment(self):
        str1 = ""
        str2 = ""
        i = len(self.backtrack_indices) - 1
        j = len(self.backtrack_indices[0]) - 1
        while i > 0:
            while j > 0:
                if (i - 1 == self.backtrack_indices[i][j][0] and
                        j - 1 == self.backtrack_indices[i][j][1]):  # diagonal
                    str1 += self.seq1[i - 1] + ''
                    str2 += self.seq2[j - 1] + ''
                elif i - 1 == self.backtrack_indices[i][j][0]:  # up direction
                    str1 += self.seq1[i - 1] + ''
                    str2 += '-'
                else:
                    str1 += '-'
                    str2 += self.seq2[j - 1] + ''  # left direciton
                k = i
                i = self.backtrack_indices[i][j][0]
                j = self.backtrack_indices[k][j][1]
        str1 = str1[::-1]
        str2 = str2[::-1]
        # Print alignment and score
        print("\n********** Result **********\n")
        print("Sequence 1: ", str1[1:])
        print("Sequence 2: ", str2[1:])
        print("Score     : ", self.middle_matrix[-1][-1])


"""MAIN"""
# Get file name from user
while True:
    input_file = input('Enter the name of the input File: ')
    if not os.path.isfile(input_file):
        print('File not found. Try again.')
        continue
    break
file = open(input_file, 'r')
# format dna as lines and make all uppercase
sequences = [line.strip() for line in file]
sequences = [seq.upper() for seq in sequences]
seq1 = sequences[0]
seq2 = sequences[1]
seq1 = " " + seq1
seq2 = " " + seq2

sequenceAlignment = SequenceAlignment(seq1, seq2)
sequenceAlignment.runAlignment()
