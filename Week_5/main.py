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


# Counting bloom filter

def addElementCounting(elt):
    for seed in seeds:
        index = mmh3.hash(str(elt), seed, signed = False) %sizeMax
        bitsArr[index] += 1

def containsElementCounting(elt):
    for seed in seeds:
        index = mmh3.hash(str(elt), seed, signed = False) % sizeMax
        if bitsArr[index] == 0:
            return False
    return True

def removeCounting(elt):
    if containsElementCounting(elt):
        for seed in seeds:
            index = mmh3.hash(str(elt), seed, signed = False) % sizeMax
            bitsArr[index] = max(0, bitsArr[index] - 1)





def main():
    addElement("Hello")
    addElement("GoodBye")
    print(containsElement("Hello"))
    print(containsElement("GoodBye"))
    print(containsElement("bonjour"))


    addElementCounting("hello")
    addElementCounting("hi")
    addElementCounting("hello")   

    print(containsElementCounting("hello"))  # True
    print(containsElementCounting("hi"))   # True
    print(containsElementCounting("bonjour")) # False

    removeCounting("hello")
    print(containsElementCounting("hello"))  # True

    removeCounting("hello")
    print(containsElementCounting("hello"))  # False



if __name__ == "__main__":
    main()
