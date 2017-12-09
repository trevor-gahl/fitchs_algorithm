import sys


sequence = []

filename = "phylogenyFile.txt"

with open("phylogenyFile.txt", 'r', encoding='utf-8') as newFile:  # open file with extended ascii
    data = newFile.readlines()
    for line in data:
        line.rstrip('\n')
        line = line.replace(' ', " ")     # strip out any whitespace
        sequence.append(line)

for x in range(len(sequence)):
    print(sequence[x])
