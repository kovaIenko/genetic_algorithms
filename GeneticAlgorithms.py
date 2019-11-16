from collections import Counter

l = 8
N = 100

# RWS
def selRoulette(list):
    sum = sumHealth(list)
    probability = []
    for i, ch in enumerate(list):
        prob = healthOf(ch)/sum
        probability.insert(i, prob)
  # we have the list of probabilities

    return list

### execute the mutation
### indexes = list of pointers to chromosomes which were mutated
def mutation(list, percentage):
    count = list.size*l
    indexes = []
    #generate the indexes which chromosomes are mutated
    #i think we don't need to save the indexes of gens which we changed
    return (list, indexes)


# generate the population by the 1 method
def init1(lenCh, numbs):
    print()

# generate the population by 2 method
def init2(lenCh, numbs):
    print()

# generate the population by the 3 method
def init3(lenCh, numbs):
    print()


# calculate the sum of all population
def sumHealth(list):
    sum = 0
    for ch in list:
        sum += healthOf(ch)
    return sum

def healthOf(ch):
   return ch.count('0')

# we do it on the beginning and after mutation
def calcHealth(list, iteration, listMutated):

    # if only we start a lifecycle
    if iteration == 0:
        for i in list:
            print()

# calculate the number different unique distances [ 1:3, 2:7, ...]
def calcNumbDistances(list):
    map = Counter(list)
    return dict(map)

# distance between two chromosome
# probably it can be changed, method compares the id of characters of chromosomes
def calcDistance(a, b):
    dim = 0
    for x, y in zip(a, b):
            if id(x) != id(y):
                dim += 1
    return dim

def execution(iterations):
    list = init1(l, N)
    print(iterations)


# 0 -> 1 or 1 -> 0
def turnOverGen(gen):
    return 0 if gen == 1 else 1

execution(100)





