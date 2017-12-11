#########################################################################
## Name: Trevor Gahl, David Schwer, Procassius Loftin                  ##
## Course: Computational Biology                                       ##
## Date: 9/25/17                                                       ##
## About: Python implementation of phylogenetic prediction             ##
#########################################################################

import sys
import binaryModule
import statistics
from binaryModule import Node, pprint, inspect, _bst_insert, _new_node, _validate_tree, _build_tree, _weight_of, _add_left, _add_right, _build_list, _right_of, _left_of, _value_of, _null
match = 2
other = -1
maxScore = 0
maxPosition = (0, 0)
pairwise_alignment_score = 0
sequence = []
finalSequence = []
score_sequence = []
deltaMatrix = []
fitchIndex = []
filename = "phylogenyFile.txt"

########################################
## Main method to run local alignment ##
########################################


def main():
    # Sequence input structure
    sequence, distance_matrix = fileReader(filename)
    for i in range(len(sequence)):
        finalSequence.append(sequence[i].strip('\n'))
    # Initialize distance_matrix
    #output = compare(finalsequence[1], sequence[2])
    # print(output)
    distance_matrix = pairwiseDistanceMatrix(finalSequence, distance_matrix)
    score_sequence, min_sequence = scoreSequence(distance_matrix)
    #print(score_sequence, min_sequence)
    # print_phylogenyTree(score_sequence)
    #printIndexAndScore(sequence, score_sequence)
    new_root = buildTree(score_sequence, finalSequence)
    pprint(new_root)
    fitchsIndexCreation(new_root)
    pprint(new_root)
    # print_phylogenyTree(deltaMatrix)

#########################################################
## Prints out the sequence and the corresponding score ##
#########################################################


def fileReader(filename):
    # open file with extended ascii
    with open(filename, 'r', encoding='utf-8') as newFile:
        data = newFile.readlines()
        for line in data:
            line.rstrip('\n')
            line = line.replace(' ', " ")     # strip out any whitespace
            sequence.append(line)
    distance_matrix = [
        [0 for col in range(len(sequence))]for row in range(len(sequence))]
    return sequence, distance_matrix


def delta(leftSequence, rightSequence):
    deltaScore = 0

    if(leftSequence == rightSequence):
        deltaScore = 0
    else:
        deltaScore = 1
    print(deltaScore)
    # print("Deltas")
    deltaMatrix.append(deltaScore)
    '''for x in range(len(deltaMatrix)):
        print(deltaMatrix[x])
        '''


# def printIndexAndScore(sequence, score_sequence):
    # for x in range(len(sequence)):
    # print(score_sequence[x])
    # print(sequence[x])
    # if(x + 1 < len(sequence)):
    #delta(sequence[x], sequence[x + 1])
    # else:
    #print("End of Sequence")


def buildTree(weights, values):
    medianValues = list(weights)
    # print(weights)
    medianValues.sort()
    # print(weights)
    # print(medianValues)
    if len(medianValues) % 2 == 0:
        root_value = medianValues[int(len(medianValues) / 2)]
    else:
        root_value = statistics.median(medianValues)

    # print(root_value)
    # print(weights)
    for x in range(len(weights)):
        if weights[x] == root_value:
            root_index = x
        else:
            continue
    new_root = _new_node(values[root_index], weights[root_index])
    for x in range(len(weights)):
        if x != root_index:
            _bst_insert(new_root, values[x], weights[x])
        else:
            continue
    return(new_root)


def fitchsIndexCreation(root_node):
    global fitchIndex
    node = root_node
    if root_node == _null:
        return
    if _left_of(node) == _null and _right_of(node) == _null:
        return
    try:
        test_1 = set(_value_of(_left_of(node)))
    except:
        test_1 = set()
    try:
        test_2 = set(_value_of(_right_of(node)))
    except:
        test_2 = set()
    inSet = set.intersection(test_1, test_2)
    if len(inSet) == 0:
        node.value = set.union(test_1, test_2)
    else:
        node.value = inSet

    fitchsIndexCreation(_left_of(node))
    fitchsIndexCreation(_right_of(node))


def compare(a, b):
    notEqual = []
    for x, y in zip(a, b):
        if x == y:
            print('equal')
            notEqual.append((x, y))
        else:
            print('not equal')
            notEqual.append((x, y))
    if len(notEqual) < 1:
        return a
    else:
        return notEqual

        ############################################
        ## Calculates the pairwise distance score ##
        ## for all sequences in S                 ##
        ############################################


def pairwiseDistanceMatrix(sequence_list, distance_matrix):
    global seq1
    global seq2
    global rows
    global cols
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
            # print(rows, cols)
            score_matrix, start_pos = createScoreMatrix(rows, cols)
            distance_score = traceback(score_matrix, start_pos)
            # print(distance_score)
            distance_matrix[x][y] = distance_score
        # offset += 1
    return distance_matrix


##############################################
## Creates the sequence of alignment scores ##
##############################################

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
    maxPosition = (0, 0)
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


######################################
## Pretty prints the phylogeny tree ##
######################################
'''
def makeTree(matrixTree,sequence0):

    for x in range(len(matrixTree)):
        ourTree.addNode(ourTree.root, matrixTree[x], sequence0[x])
        print("score and sequence")
        print(matrixTree[x])
        print(sequence0[x])
        print("end score and sequence")
        print("our tree type ")
        print(type(ourTree.root) )
        print( "test tree type")
        print(type(testTree.root))

    bt.printInorder(ourTree.root)

        makeTree(score_sequence,sequence)
        '''


if __name__ == '__main__':
    sys.exit(main())
