import mmh3
import math
import random 

sizeMax = 1000000
nbHashes = 7

bitsArr = [0] * sizeMax
seeds =[]

for _ in range(nbHashes):
    seeds.append(random.randint(0, 1000))

def addElement(elt):
    for seed in seeds:
        index = mmh3.hash(str(elt), seed, signed=False) % sizeMax
        bitsArr[index] = 1

def containsElement(elt):
    for seed in seeds:
        index = mmh3.hash(str(elt), seed, signed=False) % sizeMax
        if bitsArr[index] == 0:
            return False
    return True

addElement("Hello")
addElement("GoodBye")

print(containsElement("Hello"))
print(containsElement("GoodBye"))
print(containsElement("bonjour"))
