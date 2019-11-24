import random

import numpy as np  # linear algebra
from collections import Counter
from decimal import Decimal
from operator import attrgetter


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
    list = sorted(list, key=healthOf)
   # print(list)
    probability = []
    chosen = []
    for ch in list:
        prob = healthOf(ch) / sumH
        probability.append(prob)
  #  print(probability)
    cumsum = np.cumsum(probability)
  #  print(cumsum)
    for i in list:
        u = random.random()
     #   print("u: "+ str(u))
        for idx, val in enumerate(cumsum):
            if val >= u:
          #      print(idx)
          #      print(list[idx])
                chosen.append(list[idx].copy())
                break
    #print(chosen)
    return chosen

'''
l = 4
N = 4
pop = generate_population(l, N, 0, 1)
print(selRoulette(pop))
print("----------------")
'''


# 0 -> 1 or 1 -> 0
def turnOverGen(ch, indOfGen):
    if ch[indOfGen] == 1:
        ch[indOfGen] = 0
    else:
        ch[indOfGen] = 1
    return ch


def mutation(list, percentage, l):
    #print(percentage)
    # pm = px + 0.2*px
    # pm = px - 0.2*px
    # pm = px/2
    # pm = px/10
    # pm = px/100
    count_ = len(list) * l
    #print(count_)
    numbOfMut = count_ * percentage
    #print("NumMut: ")
    #print(numbOfMut)
    indexes = []
    # generate the indexes which chromosomes will be mutated
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

def hamming(chaine1, chaine2):
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))


# return the map [ distance : frequency ]
def calc_all_distances(list, N, l):
    frequency = { new_list: 0 for new_list in range(l + 1)}
    for i in range(0, N):
        for j in range(i + 1, N):
            dis = hamming(list[i], list[j])
            frequency[dis] = frequency.get(dis) + 1
    return frequency


# print(calc_all_distances(["01010", "10101"], 2, 5))


def healthMean(list, N):
    return  sum_health(list) / N


### execute the mutation
### indexes = list of pointers to chromosomes which were mutated

def calculate_individual_fitness(chromosome):
    return sum(map(lambda x: x == 0, chromosome))


# We select random t chromosomes, compare them with each other, choose the fittest one and send it to the mating pool. Then chromosomes are returned to the initial population
# The process continues until we select N chromosomes for a new population
def tournament_selection(population, t):
    mating_pool = []
    N = len(population)
    for chromosome in population:
        random_chromosomes = []
        health_of_random = []
        for i in range(0, t):
            rand_index = np.random.randint(0, N)
            random_chromosomes.append(population[rand_index])
            health_of_random.append(calculate_individual_fitness(random_chromosomes[i]))
        ##print(random_chromosomes)
        index_of_fittest = health_of_random.index(max(health_of_random))
        mating_pool.append(random_chromosomes[index_of_fittest])
    return mating_pool

# Testing:
##print(tournament_selection([[0, 1, 1], [1, 1, 1]], 2))

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


def save_to_file(worksheet, dict, iterate):
    print(iterate)
    row = int(iterate / CONST_NUMB_GETTING_INFO)
    col = 0
    if row == 0:
      for key in dict.keys():
         worksheet.write(row, col, key)
         col = col + 1
      row += 1
    col = 0
    for key in dict.keys():
        worksheet.write(row, col, dict.get(key))
        col = col + 1




CONST_STOP_ALGORITHM = 200000
CONST_STOP_ALGORITHM_BY_MEAN_HEALTH = 0.0001
CONST_NUMB_GETTING_INFO = 10


def should_be_stopped(worksheet, pop, N, l, i, sum_mean_health):
    if i + 1 > CONST_STOP_ALGORITHM:
        print("Algorithm was stopped")
        return True
    if i % CONST_NUMB_GETTING_INFO == 0 and i != 0:  # each 11n iterations
        save_to_file(worksheet, calc_all_distances(pop, N, l), i)
        current_mean_health = round(Decimal(healthMean(pop, N)), 4)
        print("curr")
        print(current_mean_health)
        print("last")
        last_10_mean_health = round(sum_mean_health / CONST_NUMB_GETTING_INFO, 4)
        print(last_10_mean_health)
        if current_mean_health - last_10_mean_health > CONST_STOP_ALGORITHM_BY_MEAN_HEALTH:
            print("Algorithm was stopped")
            return True
    return False





def execution(l, N):
    workbook = xlsxwriter.Workbook('data.xls')
    worksheet = workbook.add_worksheet()

    pop = generate_population(l, N, 0, 1)
    print(pop)
    pm = 1 / (10 * l)
    sum_mean_health = round(Decimal(0), 4)   #healthMean(pop, N)
    for i in range(N):
        # =================================================================================
         #print("true " + str(i) + "   " + str(np.round(Decimal(healthMean(pop, N)), 4)))

        if should_be_stopped(worksheet, pop, N, l, i, sum_mean_health):
            sum_mean_health = 0
            break

        sum_mean_health = sum_mean_health + round(Decimal(healthMean(pop, N)), 4)

        # =================================================================================
        pop = selRoulette(pop)
        print("after roulette")
        print(pop)
        pop, chs_were_muted = mutation(pop, pm, l)
        print("after mutation")
        print(chs_were_muted)
        print(pop)
        print("- - -  - - - - - - - -- - - - - - - -- - -- -- - - - - - - - ")
    workbook.close()


l = 10
N = 21
execution(l, N)
