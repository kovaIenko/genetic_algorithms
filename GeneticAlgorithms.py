import random

import numpy as np  # linear algebra
from collections import Counter
from decimal import *
from operator import attrgetter
getcontext().prec = 6


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


def healthOf(ch):
    return np.count_nonzero(ch == 0)

# RWS

def selRoulette(list, list_of_health):
    sumH = sum(list_of_health)
    list = sorted(list, key=healthOf)
    probability = []
    chosen = []
    for ch in list:
        prob = Decimal(healthOf(ch))/Decimal(sumH)
        probability.append(prob)
    cumsum = np.cumsum(probability)
    cumsum[len(cumsum)-1] = Decimal(1.)
    for i, val in enumerate(list):
        u = random.random()
        for idx, val in enumerate(cumsum):
            if val >= u:
                list_of_health[i] = healthOf(list[idx])
                chosen.append(list[idx].copy())
                break
    return chosen, list_of_health

# 0 -> 1 or 1 -> 0
def turnOverGen(ch, indOfGen):
    if ch[indOfGen] == 1:
        ch[indOfGen] = 0
    else:
        ch[indOfGen] = 1
    return ch


def mutation(list, percentage, l, list_of_health):
    count_ = len(list) * l
    numbOfMut = count_ * percentage
    indexes = []
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
        currentCh = turnOverGen(currCh, indOfGen)
        list[indOfCh] = currentCh
        list_of_health[indOfCh] = healthOf(currentCh)
        indexes.append(indOfCh)
    return list, list_of_health


# generate the population by the 3 method
def init3(lenCh, numbs):
    print()

def hamming(chaine1, chaine2):
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))


# print(calc_all_distances(["01010", "10101"], 2, 5))


def healthMean(N, list_health):
    return Decimal(sum(list_health)/N)


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


CONST_STOP_ALGORITHM = 2000
CONST_STOP_ALGORITHM_BY_MEAN_HEALTH = 0.0001
CONST_NUMB_GETTING_INFO = 10


def should_be_stopped(worksheet, pop, N, l, i, sum_mean_health, current_mean_health):
    if i % CONST_NUMB_GETTING_INFO == 0 and i != 0:  # each 11n iterations
        save_to_file(worksheet, calc_all_distances(pop, N, l), i)
        last_10_mean_health = sum_mean_health / CONST_NUMB_GETTING_INFO
        if Decimal(current_mean_health) - Decimal(last_10_mean_health) < CONST_STOP_ALGORITHM_BY_MEAN_HEALTH:
            print("the_best_individual_distances")  # wild type
            the_best_individual_distances = distances_for_wild_type(pop)
            print(the_best_individual_distances)
            print("Algorithm was stopped: there is mistake ")
            return True
    if i+1 == CONST_STOP_ALGORITHM:
        print("the_best_individual_distances")  # wild type
        the_best_individual_distances = distances_for_wild_type(pop)
        print(the_best_individual_distances)
    return False


def init_health_list(list):
    list_of_health = []
    for ch in list:
        list_of_health.append(healthOf(ch))
    return list_of_health

# return the map [ distance : frequency ]
def calc_all_distances(list, N, l):
    frequency = { new_list: 0 for new_list in range(l + 1)}
    for i in range(0, N):
        for j in range(i + 1, N):
            dis = hamming(list[i], list[j])
            frequency[dis] = frequency.get(dis) + 1
    return frequency


def calc_hamming_to_ideal(list_health, l):
    frequency = {new_list: 0 for new_list in range(l + 1)}
    for ch_health in list_health:
        dis = l-ch_health
        frequency[dis] = frequency.get(dis) + 1
        return frequency


def distances_for_wild_type(last_pop):
     frequencies = {new_list: 0 for new_list in range(l + 1)}
     frequency = [0]*l
     for ind in last_pop:
        for i, gen in enumerate(ind):
          if gen == 0:
            frequency[i] += 1
          else:
            frequency[i] -= 1
     the_best= []
     for f in frequency:
         if f >= 0:
             the_best.append(0)
         else:
             the_best.append(1)
     for ch in last_pop:
         dis = hamming(ch, the_best)
         frequencies[dis] = frequencies.get(dis) + 1
     return frequencies


def bestHealth(list_of_health):
    return Decimal(max(list_of_health))

def deviation_meanHealth_and_optimum(list_of_health, N, l ):
    mean = healthMean(N, list_of_health)
    return Decimal(l) - mean

def deviation_bestHealth_and_optimum(list_of_health, l):
    return bestHealth(list_of_health) - l


import xlsxwriter


def save_to_file(worksheet, dict, iterate):
    row = int(iterate / CONST_NUMB_GETTING_INFO) - 1
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

def save_to_file_each_iter(worksheet, iteration, list_of_health, N, l):
    meanHealth = healthMean(N, list_of_health)
    Healthbest = bestHealth(list_of_health)
    val_dev_meanHealth = deviation_meanHealth_and_optimum(list_of_health, N, l)
    val_dev_bestHealth = deviation_bestHealth_and_optimum(list_of_health, l)
    iteration




def execution(l, N):
    workbook = xlsxwriter.Workbook('data.xls')
    worksheet = workbook.add_worksheet()
    pop = generate_population(l, N, 0, 1)
    print(distances_for_wild_type(pop))
    pm = 1 / (10 * l)
    sum_mean_health = Decimal(0)
    list_health = init_health_list(pop)
    for i in range(CONST_STOP_ALGORITHM):
        mean = healthMean(N, list_health)
        if should_be_stopped(worksheet, pop, N, l, i, sum_mean_health, mean):
            sum_mean_health = 0
            break
        sum_mean_health += mean
        pop, list_health = selRoulette(pop, list_health)
        print("after roulette")
        print(pop)
        pop, list_health = mutation(pop, pm, l, list_health)
        print("after mutation")
        print(pop)
        print("- - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
    workbook.close()


l = 10
N = 300
execution(l, N)


def execution2(l, N, X ,Y):
    workbook = xlsxwriter.Workbook('data.xls')
    worksheet = workbook.add_worksheet()
    pop = generate_population(l, N, X, Y)
    pm = 1 / (10 * l)
    sum_mean_health = Decimal(0)
    list_health = init_health_list(pop)
    for i in range(CONST_STOP_ALGORITHM):
        mean = healthMean(N, list_health)
        if should_be_stopped(worksheet, pop, N, l, i, sum_mean_health, mean):
            sum_mean_health = 0
            break
        sum_mean_health += mean
        pop, list_health = selRoulette(pop, list_health)
        print("after roulette")
        print(pop)
        pop, list_health = mutation(pop, pm, l, list_health)
        print("after mutation")
        print(pop)
        print("- - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
    workbook.close()

l = 10
N = 2000
X = 0.9
Y = 0.1
#execution2(l, N, X, Y)