import random

import numpy as np  # linear algebra
from collections import Counter
from decimal import *
from operator import attrgetter

getcontext().prec = 4
import scipy.spatial.distance
import timeit  # for measuring time
import matplotlib.pyplot as plt


class UnknownCoding(Exception):
    pass


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
        raise UnknownCoding(
            "The generation function hasn't been implemented for '{0}' chromosomes coding".format(coding))
    return general_population


# Calculate the health of all population (alternative method)
def population_health(pop):
    return sum(map(calculate_individual_fitness, pop))


def healthOf(ch):
    return np.count_nonzero(ch == 0)


# RWS

def selRoulette(list, list_of_health):
    sumH = sum(list_of_health)
    list = sorted(list, key=healthOf)
    probability = []
    chosen = []
    for ch in list:
        prob = Decimal(healthOf(ch)) / Decimal(sumH)
        probability.append(prob)
    cumsum = np.cumsum(probability)
    cumsum[len(cumsum) - 1] = Decimal(1.)
    for i, val in enumerate(list):
        u = random.random()
        for idx, val in enumerate(cumsum):
            if val >= u:
                list_of_health[i] = healthOf(list[idx])
                chosen.append(list[idx].copy())
                break
    return chosen, list_of_health


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


def mutation(list, percentage, l, list_of_health):
    # print(percentage)
    # pm = px + 0.2*px
    # pm = px - 0.2*px
    # pm = px/2
    # pm = px/10
    # pm = px/100
    count_ = len(list) * l
    # print(count_)
    numbOfMut = count_ * percentage
    # print("NumMut: ")
    # print(numbOfMut)
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
        currentCh = turnOverGen(currCh, indOfGen)
        list[indOfCh] = currentCh
        list_of_health[indOfCh] = healthOf(currentCh)
        indexes.append(indOfCh)
    return list, list_of_health


# generate the population by the 3 method
def init3(lenCh, numbs):
    print()


'''
Comparing two binary strings of equal length, Hamming distance is the number of bit positions in which the two bits are different.
In order to calculate the Hamming distance between two strings, we perform their XOR operation and then count the total number of 1s in the resultant string.
'''


def hamming(chaine1, chaine2):
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))


# Alternative method (using a third-party module)
def hamming_distance(a, b):
    if len(a) != len(b):  # Checking if chromosomes have equal length
        return None
    overall_genes = len(a)
    ratio_of_different_genes = scipy.spatial.distance.hamming(a, b)
    number_of_different_genes = int(ratio_of_different_genes * overall_genes)
    return number_of_different_genes


# Return the dictionary of the form { distance : frequency }
def calc_all_distances(list, N, l):
    frequency = {new_list: 0 for new_list in range(l + 1)}
    for i in range(0, N):
        for j in range(i + 1, N):
            dis = hamming(list[i], list[j])
            frequency[dis] = frequency.get(dis) + 1
    return frequency


# print(calc_all_distances(["01010", "10101"], 2, 5))

# Mean health of the population

def healthMean(N, list_health):
    return sum(list_health) / N


# Alternative method

def mean_health_of_population(pop):
    return population_health(pop) / len(pop)


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


# Save a histogram to png file with all the parameters specified (except mutation)
def build_histogram(pop, N, l, iter_num, x, y):
    distances = calc_all_distances(pop, N, l)
    plt.bar(list(distances.keys()), distances.values(), color='g', width=0.9)
    plt.xticks(list(distances.keys()))
    plt.savefig("N={0}_l={1}_iter={2}_X={3}%_Y={4}%.png".format(N, l, iter_num, x, y))


# Build a line plot with mean health for each iteration and save to png

def build_line_graph(health_values, N, l, x, y):
    plt.clf()
    plt.plot(health_values)
    plt.savefig("health-over-generations_N={0}_l={1}_X={2}%_Y={3}%.png".format(N, l, x, y))


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


CONST_STOP_ALGORITHM = 200000
PRECISION = 0.0001
CONST_NUMB_GETTING_INFO = 10


def should_be_stopped(worksheet, pop, N, l, i, sum_mean_health, current_mean_health):
    if i % CONST_NUMB_GETTING_INFO == 0 and i != 0:  # each 11n iterations
        # save_to_file(worksheet, calc_all_distances(pop, N, l), i)
        last_10_mean_health = sum_mean_health / CONST_NUMB_GETTING_INFO
        if Decimal(current_mean_health) - Decimal(last_10_mean_health) < PRECISION:
            print("Algorithm was stopped: there is mistake ")
            return True
    return False


def init_health_list(list):
    list_of_health = []
    for ch in list:
        list_of_health.append(healthOf(ch))
    return list_of_health


def execution(l, N):
    workbook = xlsxwriter.Workbook('data.xls')
    worksheet = workbook.add_worksheet()
    pop = generate_population(l, N, 0, 1)
    pm = 1 / (10 * l)

    sum_mean_health = Decimal(0)
    list_health = init_health_list(pop)

    for i in range(CONST_STOP_ALGORITHM):
        mean = Decimal(healthMean(N, list_health))
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
        print("- - -  - - - - - - - -- - - - - - - - - - - - - - - - - - - - - ")
    workbook.close()


l = 10
N = 200000
execution(l, N)


def run_genetic_algorithm(l, N, X, Y):
    pop = generate_population(l, N, X, Y)
    counter = 0
    x = X * 100
    y = Y * 100
    health_during_generations = []
    for i in range(1, CONST_STOP_ALGORITHM + 1):
        health_of_each_individual = list(map(calculate_individual_fitness, pop))
        print("# Iteration number {0}".format(i))
        # print("Current population:")
        # print(pop)
        # print("Health of each chromosome: ")
        # print(health_of_each_individual)
        previous_mean_health = mean_health_of_population(pop)
        print("Previous mean health: ")
        print(previous_mean_health)
        #if i % 10 == 0:
            #build_histogram(pop, N, l, i, x, y)
        pop = selRoulette(pop)
        pop, mutated_chromosomes = mutation(pop, 0.01, l)
        # print("Mutated: ")
        # print(mutated_chromosomes)
        # print("Unique values: ")
        chromosomes_numbers = set(mutated_chromosomes)
        # print(set(mutated_chromosomes))
        print("{0} of {1} chromosomes has mutated".format(len(chromosomes_numbers), N))
        pop = np.asarray(pop)  # convert list to numpy array
        # Recalculate health of mutated chromosomes only
        for j in chromosomes_numbers:
            health_of_each_individual[j] = calculate_individual_fitness(pop[j])
        current_mean_health = np.mean(health_of_each_individual)
        print("Mean health after mutation:")
        print(current_mean_health)
        health_during_generations.append(current_mean_health)
        if np.abs(previous_mean_health - current_mean_health) < PRECISION and i > 1:
            counter = counter + 1
        else:
            counter = 0
        if counter >= 10:
            print("Mean health stays the same during 10 iterations -> algorithm was stopped")
            print("Iteration number: ", counter)
            build_line_graph(health_during_generations, N, l, x, y)
            break
        print("Go to the next iteration")
        if i == CONST_STOP_ALGORITHM:
            build_line_graph(health_during_generations, N, l, x, y)
        # print()


#run_genetic_algorithm(10, 100, 1, 0)