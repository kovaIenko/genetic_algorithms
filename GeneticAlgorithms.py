import random
import Properties
import xlsxwriter
import numpy as np  # linear algebra
from collections import Counter
from decimal import *


getcontext().prec = 6

import scipy.spatial.distance
import matplotlib.pyplot as plt
from math import floor

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


'''
# Calculate the health of all population (alternative method)
def population_health(pop):
    return sum(map(calculate_individual_fitness, pop))
'''


def health_of_1(ch):
    return np.count_nonzero(ch == 0)

# the health func for third type of evaluation
def health_of_3(type_, l, k):
    return {
         type_ == Properties.CONST_TYPE_NEUTRAL: l,
         type_ == Properties.CONST_TYPE_PATHOGENIC: l/(l-k),
         type_ == Properties.CONST_TYPE_LETHAL: 0.1,
    }.get(type_)

def health_of_2(l):
    return l

# RWS
def selRoulette(pop, health_func, list_of_health, N):
    chosen = []
    pop = sorted(pop, key=health_func)
    probability, list_of_health = calc_probabilities(list_of_health, N)
    cumsum = np.cumsum(probability)
    cumsum[len(cumsum) - 1] = Decimal(1.)
    list_health_temp = []
    for i in range(len(pop)):
        u = random.random()
        for idx, val in enumerate(cumsum):
            if val >= u:
                if list_of_health:
                   list_health_temp.append(list_of_health[idx])
                chosen.append(pop[idx].copy())
                break
    return chosen, list_health_temp


def calc_probabilities(health_list, N):
    probability = []
    if health_list:
        sum_health = sum(health_list)
        health_list = sorted(health_list)
        for health in health_list:
            prob = Decimal(health) / Decimal(sum_health)
            probability.append(prob)
    else:
        probability = [1 / N] * N
    return probability, health_list

# 0 -> 1 or 1 -> 0
def turnOverGen(ch, indOfGen):
    if ch[indOfGen] == 1:
        ch[indOfGen] = 0
    else:
        ch[indOfGen] = 1
    return ch


def mutation(pop, percentage, l, list_of_health, health_of, features=None, pathogenic_muted=None):
    count_ = len(pop) * l
    numbOfMut = count_ * percentage
    for i in range(int(numbOfMut)):
        index_of_ch, index_of_gen = generate_ch_and_gen(count_)
        if not features:
             currCh = pop[index_of_ch]
             currentCh = turnOverGen(currCh, index_of_gen)
             pop[index_of_ch] = currentCh
             list_of_health[index_of_ch] = health_of(currentCh)
        elif not list_of_health:
             currCh = pop[index_of_ch]
             currentCh = turnOverGen(currCh, index_of_gen)
             pop[index_of_ch] = currentCh
        else:
             current_map = features[index_of_ch]
             type = getType(current_map, index_of_gen)
             pathogenic_muted = increment_pathogenic_counter(pathogenic_muted, type)
             list_of_health[index_of_ch] = health_of(type, l)
    return pop, list_of_health, pathogenic_muted


def increment_pathogenic_counter(counter, type):
    if type == Properties.CONST_TYPE_PATHOGENIC:
        counter += 1
    return counter

def generate_ch_and_gen(count_):
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
    return indOfCh, indOfGen

def getType(current_map, ind_of_gen):
     for (type, list_of_genes) in current_map.keys():
         if ind_of_gen in list_of_genes:
             return type
     return None

#list of maps for each chromosome
def list_features_of_ch(N, l):
    maps = []
    for ch in len(N):
       maps.append(mutation_genes_distribution(l))
       maps.append(None)
    return maps

# Returns a map for a chromosome, where key is a mutation type and value is a list of genes positions
def mutation_genes_distribution(chr, first_neutral_percent=0.135, any_neutral_percent=0.245,
                                pathogenic_percent=0.0232):
    l = len(chr)
    # Calculate number of mutations of each type
    first_neutral_mutations = floor(first_neutral_percent * l)
    any_neutral_mutations = floor(any_neutral_percent * l)
    pathogenic_mutations = floor(pathogenic_percent * l)
    lethal_mutations = l - (first_neutral_mutations + any_neutral_mutations + pathogenic_mutations)
    overall_neutral = first_neutral_mutations + any_neutral_mutations
    overall = overall_neutral + pathogenic_mutations + lethal_mutations
    # Check the number of mutations of each type
    general_scheme = {'neutral_first': first_neutral_mutations, 'neutral_other': any_neutral_mutations, 'pathogenic':
        pathogenic_mutations, 'lethal': lethal_mutations}

    specific_scheme = dict()
    first_neutral_indices = list(range(first_neutral_mutations))
    print("First neutral indices:")
    print(first_neutral_indices)


    other_neutral_indices = []
    # Locate neutral locuses

    i = first_neutral_mutations # choose randomly 24.5 % of other neutral genes (if l = 100, i is in range [13, 37) = > we get 24 genes)
    # We need to iterate until we locate all neutral locuses
    while i != overall_neutral:
        rand = np.random.randint(first_neutral_mutations, l)  # if l = 100, a random number from [13, 100)
        if rand in other_neutral_indices:  # if there's a collision, then go one iteration back
            continue
        else:
            other_neutral_indices.append(rand)
            i = i + 1
    print("Other neutral indices:")
    print(other_neutral_indices)
    neutral_indices = np.concatenate((first_neutral_indices, other_neutral_indices), axis=0)
    print("Indices of neutral locuses:")
    print(neutral_indices)
    print("Number of neutral indices:")
    print(len(neutral_indices))

    # Set indices in the chromosome's schemes
    specific_scheme['neutral'] = neutral_indices

    # Locate pathogenic locuses
    pathogenic_indices = []
    j = 0
    while j != pathogenic_mutations:
        rand = np.random.randint(first_neutral_mutations, l)
        if rand in neutral_indices or rand in pathogenic_indices:  # check if this locus is not already included as neutral or pathogenic
            continue
        else:
            pathogenic_indices.append(rand)
            j = j + 1
    print("Indices of pathogenic locuses:")
    print(pathogenic_indices)
    specific_scheme['pathogenic'] = pathogenic_indices


    # Locate lethal locuses
    lethal_indices = []

    k = 0
    while k != lethal_mutations:
        rand = np.random.randint(first_neutral_mutations, l)
        if rand in neutral_indices or rand in pathogenic_indices or rand in lethal_indices:
            continue
        else:
            lethal_indices.append(rand)
            k = k + 1

    specific_scheme['lethal'] = lethal_indices
    print("Lethal indices:", lethal_indices)
    print("Number of lethal indices:", len(lethal_indices))
    # Check if lists of indices don't contain equal elements (pairwise intersection must be an empty set)
    assert(set(neutral_indices) & set(pathogenic_indices) == set())
    assert(set(neutral_indices) & set(lethal_indices) == set())
    assert(set(pathogenic_indices) & set(lethal_indices) == set())
    #return general_scheme, overall, specific_scheme # for testing
    return specific_scheme


# Testing mutation_genes_distribution() method
chrom = np.random.randint(0, 2, 100)
scheme = mutation_genes_distribution(chrom)
print(scheme)


def mkdir_p(mypath):
    '''Creates a directory. equivalent to using mkdir -p on the command line'''

    from errno import EEXIST
    from os import makedirs,path

    try:
        makedirs(mypath)
    except OSError as exc: # Python >2.5
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise

'''def mutation_genes_distribution(l):
    initial = int(0.135*l)
    rest = l - initial
    list_1 = [CONST_TYPE_NEUTRAL]*initial
    list_2 = np.random.choice([CONST_TYPE_NEUTRAL, CONST_TYPE_PATHOGENIC, CONST_TYPE_LETHAL], size=rest, p=[0.2932, 0.0178, 0.689])
    result_list = np.concatenate((list_1, list_2), axis=0)
    return arrangement_list(result_list)

def arrangement_list(list_):
    dict = {CONST_TYPE_NEUTRAL: [], CONST_TYPE_PATHOGENIC: [], CONST_TYPE_LETHAL: []}
    for ind, val in enumerate(list_):
        dict[val].append(ind)
    return dict'''

#inds = mutation_genes_distribution(100);

#print(len(inds.get(CONST_TYPE_NEUTRAL)))
#print(len(inds.get(CONST_TYPE_PATHOGENIC)))
#print(len(inds.get(CONST_TYPE_LETHAL)))
'''
Comparing two binary strings of equal length, Hamming distance is the number of bit positions in which the two bits are different.
In order to calculate the Hamming distance between two strings, we perform their XOR operation and then count the total number of 1s in the resultant string.
'''


def hamming(chaine1, chaine2):
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))


'''# Alternative method (using a third-party module)
def hamming_distance(a, b):
    if len(a) != len(b):  # Checking if chromosomes have equal length
        return None
    overall_genes = len(a)
    ratio_of_different_genes = scipy.spatial.distance.hamming(a, b)
    number_of_different_genes = int(ratio_of_different_genes * overall_genes)
    return number_of_different_genes'''


# Return the dictionary of the form { distance : frequency }
def calc_all_distances(list_, N_, l_):
    frequency = {new_list: 0 for new_list in range(l_ + 1)}
    for i in range(0, N_):
        for j in range(i + 1, N):
            dis = hamming(list_[i], list_[j])
            frequency[dis] = frequency.get(dis) + 1
    return frequency

def healthMean(N, list_health):
    return Decimal(sum(list_health)/N)

'''# Alternative method

def mean_health_of_population(pop):
    return population_health(pop) / len(pop)

'''
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

import os
# Save a histogram to png file with all the parameters specified
# (a histogram with frequencies of pairwise Hamming distances between chromosomes)

def build_first_histogram(pop, N, l, iter_num, x, y, selection_type, pm):
    script_dir = os.path.dirname(__file__)
    results_dir = os.path.join(script_dir, 'FirstHistType; Selection={0},N={1};l={2};X={3};Y={4}/'.format(selection_type, N, l, x, y))
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    distances = calc_all_distances(pop, N, l)
    plt.bar(list(distances.keys()), distances.values(), color='g', width=0.9)
    plt.xticks(list(distances.keys()))
    plt.savefig(results_dir + "/iter={0};pm={1}.png".format(iter_num, pm))


# Save a histogram to png file with all the parameters specified
# (a histogram with frequencies of Hamming distances to the optimal chromosome from each chromosome in the population)

def build_second_histogram(list_health, N, l, iter_num, x, y, selection_type, pm):
    script_dir = os.path.dirname(__file__)
    results_dir = os.path.join(script_dir, 'SecondHistType; Selection={0},N={1};l={2};X={3};Y={4}/'.format(selection_type, N, l, x, y))
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    distances = calc_hamming_to_ideal(list_health, l)
    plt.bar(list(distances.keys()), distances.values(), color='g', width=0.9)
    plt.xticks(list(distances.keys()))
    plt.savefig(results_dir + "/iter={0};pm={1}.png".format(iter_num, pm))


# Save a histogram to png file with all the parameters specified
# (a histogram with frequencies of Hamming distances to the optimal chromosome from each chromosome in the population)

def build_third_histogram(pop, N, l, iter_num, x, y, selection_type, pm):
    script_dir = os.path.dirname(__file__)
    results_dir = os.path.join(script_dir, 'ThirdHistType; Selection={0},N={1};l={2};X={3};Y={4}/'.format(selection_type, N, l, x, y))
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    distances = distances_for_wild_type(pop)
    plt.bar(list(distances.keys()), distances.values(), color='g', width=0.9)
    plt.xticks(list(distances.keys()))
    plt.savefig(results_dir + "/iter={0};pm={1}.png".format(iter_num, pm))

# Build a line plot with mean health for each iteration and save to png

def build_line_graph(health_values, N, l, x, y):
    plt.clf()
    plt.plot(health_values)
    plt.savefig("health-over-generations_N={0}_l={1}_X={2}%_Y={3}%.png".format(N, l, x, y))
'''  sample = pop[:-3:-1]  # take two last chromosomes from population
  print("Chromosomes: \n", sample)
  proportion_of_different_genes = scipy.spatial.distance.hamming(sample[0], sample[1])
  print("{}% of different genes".format(proportion_of_different_genes * 100))
  print("The number of different genes is {0}".format(int(proportion_of_different_genes * l)))
 
 '''


def init_health_list(list, health_func):
    list_of_health = []
    for ch in list:
        list_of_health.append(health_func(ch))
    return list_of_health

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

def percent_polym_genes(list_of_health, l, N):
    suma = 0
    all_ = l * N
    for ch_health in list_of_health:
        suma += l - ch_health
    return suma*100/all_

def save_to_file(worksheet, dict, iterate):
    row = int(iterate / Properties.CONST_NUMB_GETTING_INFO) - 1
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

def save_to_file_the_end(worksheet, list_of_health, N, l, stop_mutated_iter = None):
    meanHealth = healthMean(N, list_of_health)
    healthbest = bestHealth(list_of_health)
    val_dev_meanHealth = deviation_meanHealth_and_optimum(list_of_health, N, l)
    val_dev_bestHealth = deviation_bestHealth_and_optimum(list_of_health, l)
    polym_genes = percent_polym_genes(list_of_health, l, N)


def should_be_stopped(worksheet, pop, N, l, i, sum_mean_health, current_mean_health):
    if i % Properties.CONST_NUMB_GETTING_INFO == 0 and i != 0:  # each 11n iterations
        # save_to_file(worksheet, calc_all_distances(pop, N, l), i)
        last_10_mean_health = sum_mean_health / Properties.CONST_NUMB_GETTING_INFO
        if Decimal(current_mean_health) - Decimal(last_10_mean_health) < Properties.PRECISION:
            print("Algorithm was stopped: there is mistake ")
            return True
    return False

# how many mutations we can do through the all iterations
def num_of_pathogenic_genes(N,l, percentage):
    return N * l * percentage

def execution(l, N):
    workbook = xlsxwriter.Workbook('data.xls')
    worksheet = workbook.add_worksheet()
    pop = generate_population(l, N, 0, 1)
    pm = 1 / (10 * l)
    sum_mean_health = Decimal(0)
    list_health = init_health_list(pop, health_of_1)
    for i in range(Properties.CONST_STOP_ALGORITHM):
        mean = healthMean(N, list_health)
        if should_be_stopped(worksheet, pop, N, l, i, sum_mean_health, mean):
            sum_mean_health = 0
            break
        sum_mean_health += mean
        pop, list_health = selRoulette(pop, health_of_1, list_health, N)
        '''print("after roulette")
        print(pop)
        print(list_health)'''
        pop, list_health, pathogenic_muted = mutation(pop, pm, l, list_health, health_of_1)
        '''print("after mutation")
        print(pop)'''
        if i % 10 == 0:
            build_first_histogram(pop, N, l, i, 0, 100, "Roulette", pm)
            build_second_histogram(list_health, N, l, i, 0, 100, "Roulette", pm)
            build_third_histogram(pop, N, l, i, 0, 100, "Roulette", pm)
        ### print("- - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
    workbook.close()


l = 8
N = 10
execution(l, N)

