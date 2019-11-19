import random

import numpy as np  # linear algebra
from collections import Counter

# generate initial population of size N. each chromosome has a length l

# X is a percentage of chromosomes that contain only 0
# Y is a percentage of chromosomes that are generated using a uniform distribution

def generate_population(l, N, X, Y, coding="binary"):
    if coding == "binary":
        x = int(X * N)
        y = int(Y * N)
        population_with_zeros = np.random.randint(0, 1, (x, l))
        uniform_population = np.random.randint(0, 2, (y, l))
        general_population = np.concatenate((population_with_zeros, uniform_population), axis=0)
    else:
        return None
    return general_population



'''
Comparing two binary strings of equal length, Hamming distance is the number of bit positions in which the two bits are different.
In order to calculate the Hamming distance between two strings, we perform their XOR operation and then count the total number of 1s in the resultant string.
'''

def compare_chromosome_fitness(first_chromosome, second_chromosome, fitness_func="hamming distance"):
    if fitness_func == "hamming distance":
        result = []
        if len(first_chromosome) == len(second_chromosome):
            length = len(first_chromosome)
            for i in range(length):
                result.append(first_chromosome[i] ^ second_chromosome[i])
            count_zeros = sum(map(lambda x: x == 0, result))
            return count_zeros
    else:
        return None

# calculate the sum of all population
def sum_health(list):
    suum = 0
    for ch in list:
        suum += healthOf(ch)
    return suum


def healthOf(ch):
    return np.count_nonzero(ch == 0)


# RWS

def selRoulette(list):
    sumH = sum_health(list)
    probability = []
    chosen = []
    for i, ch in enumerate(list):
        prob = healthOf(ch) / sumH
        probability.append(prob)

    probability.sort()

    for i in range(len(probability)):

        u = random.random()
        sum_ = 0
        for ind in probability:
            sum_ += probability[ind]
            if sum_ > u:
                chosen.append(ind)
                break

    return chosen


# 0 -> 1 or 1 -> 0
def turnOverGen(ch, indOfGen):
    chList = list(ch)
    if chList[indOfGen] == '1':
        chList[indOfGen] = '0'
    else:
        chList[indOfGen] = '1'
    return "".join(chList)


def mutation(list, percentage):
    l = 5
    px = 1 / (10 * l)
    pm = 0.1
    # pm = px + 0.2*px
    # pm = px - 0.2*px
    # pm = px/2
    # pm = px/10
    # pm = px/100
    count_ = len(list) * l
    numbOfMut = count_ * pm
    print(list)
    indexes = []
    # generate the indexes which chromosomes are mutated
    for i in range(int(numbOfMut)):
        u = int(random.random() * count_)
        if u < l - 1:
            indOfCh = 0
            if u != 0:
                indOfGen = u - 1
            else:
                indOfGen = u
        else:
            indOfCh = int((u - 1) / l)
            indOfGen = int((u - 1) % l)
        currCh = list[indOfCh]
        list[indOfCh] = turnOverGen(currCh, indOfGen)
        indexes.append(indOfCh)
    return list, indexes


# generate the population by the 3 method
def init3(lenCh, numbs):
    print()


# calculate the number different unique distances [ 1:3, 2:7, ...]
def calcNumbDistances(list):
    map = Counter(list)
    return dict(map)


def hamming(chaine1, chaine2):
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))


# return the map [ distance : frequency ]
def calc_all_distances(list, N, l):
    frequency = {new_list: 0 for new_list in range(l + 1)}
    for i in range(0, N):
        for j in range(i + 1, N):
            dis = hamming(list[i], list[j])
            frequency[dis] = frequency.get(dis) + 1
    return frequency

def healthMean(list, N):
    return sum_health(list) / N



### execute the mutation
### indexes = list of pointers to chromosomes which were mutated

def tournament_selection(t):
    if t == 2:
        pass
    elif t == 4:
        pass
    elif t == 12:
        pass


# The Hamming distance between 1-D arrays `u` and `v`, is simply the
#     proportion of disagreeing components in `u` and `v`


# test

'''  sample = pop[:-3:-1]  # take two last chromosomes from population
  print("Chromosomes: \n", sample)
  proportion_of_different_genes = scipy.spatial.distance.hamming(sample[0], sample[1])
  print("{}% of different genes".format(proportion_of_different_genes * 100))
  print("The number of different genes is {0}".format(int(proportion_of_different_genes * l)))
  '''

import xlsxwriter


def save_to_file(dict):
    print(dict)
    workbook = xlsxwriter.Workbook('data.xls')
    worksheet = workbook.add_worksheet()

    row = 0
    col = 0
    for key in dict.keys():
        worksheet.write(row, col, key)
        col = col + 1

    row = row+1
    col = 0
    for key in dict.keys():
        worksheet.write(row, col, dict.get(key))
        col = col+1

    workbook.close()


CONST_STOP_ALGORITHM = 200000
CONST_STOP_ALGORITHM_BY_MEAN_HEALTH = 0.0001


def execution(l, N):
    pop = generate_population(l, N, 0, 1)
    histogram_data_arr = []

    previous_mean_health = healthMean(pop, N)
    for i in range(N):
        current_mean_health = healthMean(pop, N)
        if i > CONST_STOP_ALGORITHM or current_mean_health - previous_mean_health > CONST_STOP_ALGORITHM_BY_MEAN_HEALTH:
            break
        #####################################
        #  EVALUATION
        #####################################
        if i > 9 and i % 10 == 0:  # there we save the data for histogram
            save_to_file(calc_all_distances(pop, N, l))


l = 3
N = 11
execution(l, N);
