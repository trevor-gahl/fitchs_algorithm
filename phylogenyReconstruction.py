import sys


sequence = [] 

filename = "phylogenyFile.txt"
#pattern = input("Enter pattern, P, that you want to find matches in string: \n")
#pattern = pattern.replace(' ',"")     # strip out any whitespace
with open("phylogenyFile.txt", 'r',encoding='utf-8') as newFile:        # if passing a file, read file and strip all whitespace and newlines
    data =newFile.readlines()
    for line in data:
        #print(data)
        line.rstrip('\n')
        sequence.append(line)
        #sequence = sequence + line.replace(' ', "")
        #print(line)
        #print("Pass no arguments for prompt to enter string and pattern, or pass one argument for file path to file containing sequence.")
#print(data)
print(sequence[0])
for x in range(len(sequence)):
    print(sequence[x])
#rint(data)
