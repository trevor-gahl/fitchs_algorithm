#########################################################################
## Name: Trevor Gahl, David Schwer, Procassius Loftin                  ##
## Course: Computational Biology                                       ##
## Date: 9/25/17                                                       ##
## About: Python implementation of phylogenetic prediction             ##
#########################################################################

import sys


match = 2
other = -1
maxScore = 0
maxPosition = (0, 0)
pairwise_alignment_score = 0
########################################
## Main method to run local alignment ##
########################################

sequence = []

filename = "phylogenyFile.txt"

# open file with extended ascii
with open("phylogenyFile.txt", 'r', encoding='utf-8') as newFile:
    data = newFile.readlines()
    for line in data:
        line.rstrip('\n')
        line = line.replace(' ', " ")     # strip out any whitespace
        sequence.append(line)

distance_matrix = [
    [0 for col in range(len(sequence))]for row in range(len(sequence))]


def main():
    # Sequence input structure

    distance_matrix = pairwiseDistanceMatrix(sequence)
    # for i in range(len(distance_matrix)):
    # print(distance_matrix[i])
    score_sequence, min_sequence = scoreSequence(distance_matrix)
    print(score_sequence, min_sequence)

    printIndexAndScore(sequence, score_sequence)


def printIndexAndScore(sequence, score_sequence):
    for x in range(len(sequence)):
        print(score_sequence[x])
        print(sequence[x])


def pairwiseDistanceMatrix(sequence_list):
    global seq1
    global seq2
    global rows
    global cols
    global distance_matrix
    for x in range(len(sequence_list)):
        for y in range(len(sequence_list)):
            # print("iteration")
            seq1 = sequence_list[x]
            seq2 = sequence_list[y]
            '''
            temp = y + offset
            if temp >= len(sequence_list):
                temp = len(sequence_list)
                # break
            else:
                seq2 = sequence_list[temp]
            '''
            rows = len(seq1) + 1
            cols = len(seq2) + 1
            #print(rows, cols)
            score_matrix, start_pos = createScoreMatrix(rows, cols)
            distance_score = traceback(score_matrix, start_pos)
            # print(distance_score)
            distance_matrix[x][y] = distance_score
        #offset += 1
    return distance_matrix


def scoreSequence(d_matrix):
    score_min = 10000000
    score_value = 0
    sc_sequence = []
    for x in range(len(d_matrix)):
        for y in range(len(d_matrix)):
            score_value = score_value + d_matrix[x][y]
        # print(score_value)
        if score_value < score_min:
            score_min_index = x
            score_min = score_value
        sc_sequence.append(score_value)
        # print(sc_sequence)
        score_value = 0
    return sc_sequence, score_min_index
###############################################
## Creates scoring matrix from input strings ##
###############################################


def createScoreMatrix(rows, cols):
    global maxScore
    score_matrix = [[0 for col in range(cols)]for row in range(rows)]

    for i in range(1, rows):
        for j in range(1, cols):
            similarity = match if seq1[i - 1] == seq2[j - 1] else other
            diag_score = score_matrix[i - 1][j - 1] + similarity
            up_score = score_matrix[i - 1][j] + other
            left_score = score_matrix[i][j - 1] + other
            curMax = max(0, diag_score, up_score, left_score)
            if curMax > maxScore:
                maxScore = curMax
                maxPosition = (i, j)

            score_matrix[i][j] = curMax
            # print(similarity)
    maxScore = 0
    return score_matrix, maxPosition

#######################################
## Creates best fit alignment string ##
## based on scoring matrix           ##
#######################################


def traceback(score_matrix, start_pos):
    pairwise_alignment_score = 0
    END, DIAG, UP, LEFT = range(4)
    aligned_seq1 = []
    aligned_seq2 = []
    x, y = start_pos
    move = nextMove(score_matrix, x, y)
    while move != END:
        if move == DIAG:
            aligned_seq1.append(seq1[x - 1])
            aligned_seq2.append(seq2[y - 1])
            x -= 1
            y -= 1
        elif move == UP:
            aligned_seq1.append(seq1[x - 1])
            aligned_seq2.append('-')
            x -= 1
            pairwise_alignment_score += 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[y - 1])
            y -= 1
            pairwise_alignment_score += 1

        move = nextMove(score_matrix, x, y)

    aligned_seq1.append(seq1[x - 1])
    aligned_seq2.append(seq1[y - 1])

    return pairwise_alignment_score

################################################
## Determines the next move for the traceback ##
################################################


def nextMove(score_matrix, x, y):

    # Assign the diagonal score
    diag = score_matrix[x - 1][y - 1]

    # Assign insertion/deletion scores
    up = score_matrix[x - 1][y]
    left = score_matrix[x][y - 1]

    # Check all three cases to find next character/insertion/deletion
    if diag >= up and diag >= left:
        return 1 if diag != 0 else 0
    elif up > diag and up >= left:
        return 2 if up != 0 else 0
    elif left > diag and left > up:
        return 3 if left != 0 else 0

    # Error detection
    else:
        # Execution should not reach here.
        raise ValueError('invalid move during traceback')

###############################################
## Creates the alignment string for printing ##
###############################################


def alignment_string(aligned_seq1, aligned_seq2):

    # Sets initial values
    idents, gaps, mismatches = 0, 0, 0
    alignment_string = []

    # Runs through both strings
    for base1, base2 in zip(aligned_seq1, aligned_seq2):

        # Checks for match
        if base1 == base2:
            alignment_string.append('|')
            idents += 1

        # Checks for insertion/deletion
        elif '-' in (base1, base2):
            alignment_string.append(' ')
            gaps += 1

        # If neither of the above, it's mismatch
        else:
            alignment_string.append(':')
            mismatches += 1

    # Returns the "alignment" string and the alignment characteristics
    return ''.join(alignment_string), idents, gaps, mismatches


if __name__ == '__main__':
    sys.exit(main())
