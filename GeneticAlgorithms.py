import random
import Properties
import csv
import numpy as np  # linear algebra
from collections import Counter
from decimal import *
import os

getcontext().prec = 6

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
        type_ == Properties.CONST_TYPE_PATHOGENIC: l / (l - k),
        type_ == Properties.CONST_TYPE_LETHAL: 0.1
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


def mutation(pop, percentage, l, list_of_health, health_of, features=None):
    count_ = len(pop) * l
    numbOfMut = count_ * percentage
    pathogenic_muted_bool = False
    for i in range(int(numbOfMut)):
        index_of_ch, index_of_gen = generate_ch_and_gen(count_, l)
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
            pathogenic_muted_bool = increment_pathogenic_counter(type)
            list_of_health[index_of_ch] = health_of(type, l)
    return pop, list_of_health, pathogenic_muted_bool


def increment_pathogenic_counter(counter, type):
    return type == Properties.CONST_TYPE_PATHOGENIC


def generate_ch_and_gen(count_, l):
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


# list of maps for each chromosome
def list_features_of_ch(N, l):
    maps = []
    for ch in len(N):
        maps.append(mutation_genes_distribution(l))
        maps.append(None)


# list of maps for each chromosome
def list_features_of_ch(pop):
    maps = []
    for ch in pop:
        # maps.append(mutation_genes_distribution(ch))
        maps.append(None)
    return maps


# Returns a map for a chromosome, where key is a mutation type and value is a list of genes positions
def mutation_genes_distribution(l, first_neutral_percent=0.135, any_neutral_percent=0.245,
                                pathogenic_percent=0.0232):
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
    # print("First neutral indices:")
    # print(first_neutral_indices)

    other_neutral_indices = []
    # Locate neutral locuses

    i = first_neutral_mutations  # choose randomly 24.5 % of other neutral genes (if l = 100, i is in range [13, 37) = > we get 24 genes)
    # We need to iterate until we locate all neutral locuses
    while i != overall_neutral:
        rand = np.random.randint(first_neutral_mutations, l)  # if l = 100, a random number from [13, 100)
        if rand in other_neutral_indices:  # if there's a collision, then go one iteration back
            continue
        else:
            other_neutral_indices.append(rand)
            i = i + 1
    # print("Other neutral indices:")
    # print(other_neutral_indices)
    neutral_indices = np.concatenate((first_neutral_indices, other_neutral_indices), axis=0)
    ''' print("Indices of neutral locuses:")
    print(neutral_indices)
    print("Number of neutral indices:")
    print(len(neutral_indices))'''

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
    # print("Indices of pathogenic locuses:")
    # print(pathogenic_indices)
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
    # print("Lethal indices:", lethal_indices)
    # print("Number of lethal indices:", len(lethal_indices))
    # Check if lists of indices don't contain equal elements (pairwise intersection must be an empty set)
    assert (set(neutral_indices) & set(pathogenic_indices) == set())
    assert (set(neutral_indices) & set(lethal_indices) == set())
    assert (set(pathogenic_indices) & set(lethal_indices) == set())
    # return general_scheme, overall, specific_scheme # for testing
    return specific_scheme


'''
# Testing mutation_genes_distribution() method
scheme = mutation_genes_distribution(100)
print(scheme)'''

'''def mutation_genes_distribution(l):
    initial = int(0.135*l)
=======
def mutation_genes_distribution(l):
    initial = int(0.135 * l)
>>>>>>> Stashed changes
    rest = l - initial
    list_1 = [CONST_TYPE_NEUTRAL] * initial
    list_2 = np.random.choice([CONST_TYPE_NEUTRAL, CONST_TYPE_PATHOGENIC, CONST_TYPE_LETHAL], size=rest,
                              p=[0.2932, 0.0178, 0.689])
    result_list = np.concatenate((list_1, list_2), axis=0)
    return arrangement_list(result_list)


def arrangement_list(list_):
    dict = {CONST_TYPE_NEUTRAL: [], CONST_TYPE_PATHOGENIC: [], CONST_TYPE_LETHAL: []}
    for ind, val in enumerate(list_):
        dict[val].append(ind)
    return dict'''

# inds = mutation_genes_distribution(100);

# print(len(inds.get(CONST_TYPE_NEUTRAL)))
# print(len(inds.get(CONST_TYPE_PATHOGENIC)))
# print(len(inds.get(CONST_TYPE_LETHAL)))
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
        for j in range(i + 1, N_):
            dis = hamming(list_[i], list_[j])
            frequency[dis] = frequency.get(dis) + 1
    return frequency


def healthMean(N, list_health):
    return Decimal(sum(list_health)) / Decimal(N)


def calculate_individual_fitness(chromosome):
    return sum(map(lambda x: x == 0, chromosome))


# We select random t chromosomes, compare them with each other, choose the fittest one and send it to the mating pool. Then chromosomes are returned to the initial population
# The process continues until we select N chromosomes for a new population
def tournament_selection(population, health_func, list_of_health, t):
    mating_pool = []
    N = len(population)
    for chromosome in population:
        random_chromosomes = []
        health_of_random = []
        for i in range(0, t):
            rand_index = np.random.randint(0, N)
            random_chromosomes.append(population[rand_index])
            health_of_random.append(health_func(random_chromosomes[i]))
            # health_of_random.append(calculate_individual_fitness(random_chromosomes[i]))
        ##print(random_chromosomes)
        index_of_fittest = health_of_random.index(max(health_of_random))
        list_of_health[rand_index] = health_of_random[index_of_fittest]
        mating_pool.append(random_chromosomes[index_of_fittest])
    return mating_pool, list_of_health


# Save a histogram to png file with all the parameters specified
# (a histogram with frequencies of pairwise Hamming distances between chromosomes)

def build_first_histogram(pop, N, l, iter_num, x, y, selection_type, pm, attempt, init_type):
    script_dir = os.path.dirname(__file__) + '/Plots'
    results_dir = os.path.join(script_dir,
                               'HammingHist;InitType={0};Attempt={1};Selection={2},N={3};l={4}/'.format(init_type + 1, attempt, selection_type, N, l))
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    distances = calc_all_distances(pop, N, l)
    plt.bar(list(distances.keys()), distances.values(), color='g', width=0.9)
    plt.xticks(list(distances.keys()))
    plt.savefig(results_dir + "/X={0};Y={1};iter={2};pm={3}.png".format(x, y, iter_num, pm))
    plt.clf()


# Save a histogram to png file with all the parameters specified
# (a histogram with frequencies of Hamming distances to the optimal chromosome from each chromosome in the population)

def build_second_histogram(list_health, N, l, iter_num, x, y, selection_type, pm, attempt, init_type):
    script_dir = os.path.dirname(__file__) + '/Plots'
    results_dir = os.path.join(script_dir,
                               'OptimalHist;InitType={0};Attempt={1};Selection={2},N={3};l={4}/'.format(init_type + 1, attempt, selection_type, N, l))
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    distances = calc_hamming_to_ideal(list_health, l)
    plt.bar(list(distances.keys()), distances.values(), color='red', width=0.9)
    plt.xticks(list(distances.keys()))
    plt.savefig(results_dir + "/X={0};Y={1};iter={2};pm={3}.png".format(x, y, iter_num, pm))
    plt.clf()


# Save a histogram to png file with all the parameters specified
# (a histogram with frequencies of Hamming distances to the optimal chromosome from each chromosome in the population)

def build_third_histogram(pop, N, l, iter_num, x, y, selection_type, pm, attempt, init_type):
    script_dir = os.path.dirname(__file__) + '/Plots'
    results_dir = os.path.join(script_dir,
                               'WildTypeHist;InitType={0};Attempt={1};Selection={2},N={3};l={4}/'.format(init_type + 1, attempt, selection_type, N, l))
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    distances = distances_for_wild_type(pop, l)
    plt.bar(list(distances.keys()), distances.values(), color='blue', width=0.9)
    plt.xticks(list(distances.keys()))
    plt.savefig(results_dir + "/X={0};Y={1};=iter={2};pm={3}.png".format(x, y, iter_num, pm))
    plt.clf()


# Build a line plot with mean health for each iteration and save to png

def build_line_graph(health_values, N, l, x, y, selection_type, pm, attempt, init_type):
    script_dir = os.path.dirname(__file__) + '/Plots'
    results_dir = os.path.join(script_dir,
                               'LineGraph;InitType={0};Attempt={1};Selection={2},N={3};l={4};/'.format(init_type + 1, attempt, selection_type, N, l))
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    plt.plot(health_values)
    plt.savefig(results_dir + "/X={0};Y={1};pm={2}.png".format(x, y, pm))
    plt.clf()


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
        dis = l - ch_health
        frequency[dis] = frequency.get(dis) + 1
    return frequency


# Ці відстані рахуємо тільки на останній ітерації (коли спрацьовує зупинка за здоров'ям, чи за кількістю ітерацій)
def distances_for_wild_type(last_pop, l):
    frequencies = {new_list: 0 for new_list in range(l + 1)}
    frequency = [0] * l
    for ind in last_pop:
        for i, gen in enumerate(ind):
            if gen == 0:
                frequency[i] += 1
            else:
                frequency[i] -= 1
    the_best = []
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


def deviation_meanHealth_and_optimum(list_of_health, N, l):
    mean = healthMean(N, list_of_health)
    return Decimal(l) - mean


def deviation_bestHealth_and_optimum(list_of_health, l):
    return bestHealth(list_of_health) - l


def percent_polym_genes(list_of_health, l, N):
    suma = 0
    all_ = l * N
    for ch_health in list_of_health:
        suma += l - ch_health
    return suma * 100 / all_


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


def save_to_file_the_end(attempt, iteration, list_of_health, N, l, pm, type_of_selection, X, Y, params_tour=-1,
                         stop_mutated_iter=None):
    row = []
    row.append(attempt)
    row.append(iteration)
    row.append(pm)
    if type_of_selection == "Tournament":
        type_of_selection = type_of_selection + ' ' + str(params_tour)
    row.append(type_of_selection)
    row.append(l)
    row.append(N)
    row.append(X)
    row.append(Y)
    row.append(healthMean(N, list_of_health))
    row.append(bestHealth(list_of_health))
    row.append(deviation_meanHealth_and_optimum(list_of_health, N, l))
    row.append(deviation_bestHealth_and_optimum(list_of_health, l))
    row.append(percent_polym_genes(list_of_health, l, N))
    if stop_mutated_iter:
        row.append(stop_mutated_iter)
    print(row)
    file = open('data.csv', 'a')
    with file:
        writer = csv.writer(file, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(row)
    file.close()


def should_be_stopped(i, sum_mean_health, current_mean_health):
    if i % Properties.CONST_NUMB_GETTING_INFO == 0 and i != 0:  # each 11n iterations
        last_10_mean_health = sum_mean_health / Properties.CONST_NUMB_GETTING_INFO
        if Decimal(current_mean_health) - Decimal(last_10_mean_health) < Properties.PRECISION:
            print("Algorithm was stopped: there is mistake ")
            return True
    return False


# how many mutations we can do through the all iterations
def num_of_pathogenic_genes(N, l, percentage):
    return N * l * percentage


'''
def execution(l, N):
    workbook = xlsxwriter.Workbook('data.numbers')
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
        print("after roulette")
        print(pop)
        print(list_health)
        pop, list_health, pathogenic_muted = mutation(pop, pm, l, list_health, health_of_1)
        print("after mutation")
        print(pop)
        if i % 10 == 0:
            build_first_histogram(pop, N, l, i, 0, 100, "Roulette", pm)
            build_second_histogram(list_health, N, l, i, 0, 100, "Roulette", pm)
            build_third_histogram(pop, N, l, i, 0, 100, "Roulette", pm)
        ### print("- - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
    workbook.close()'''

# l = 8
# N = 10
# execution(l, N)

mutation_probability = 0.001113


# At each h iterations we build only histograms for pairwise Hamming distances & Hamming distances to the optimal
# chromosome

def build_histograms(pop, N, l, i, x, y, type_of_selection, pm, health_list, attempt, init_type):
    build_first_histogram(pop, N, l, i, x, y, type_of_selection, pm, attempt, init_type)
    build_second_histogram(health_list, N, l, i, x, y, type_of_selection, pm, attempt, init_type)


'''
# For the first initialization type
def run_genetic_algorithm(l, N, X, Y, pm, type_of_selection):
    pop = generate_population(l, N, X, Y)
    counter = 0
    x = X * 100
    y = Y * 100
    health_during_generations = []
    health_list = init_health_list(pop, health_of_1)
    for i in range(1, Properties.CONST_STOP_ALGORITHM + 1):
        print("# Iteration number {0}".format(i))
        previous_mean_health = healthMean(N, health_list)
        print("Previous mean health: ")
        print(previous_mean_health)
        if i % Properties.CONST_FREQUENCY_PRINT_DIAGRAM == 0:
            build_histograms(pop, N, l, i, x, y, type_of_selection, pm, health_list)
        pop, health_list = selRoulette(pop, health_of_1, health_list, N)
        pop, health_list, pathogenic_muted = mutation(pop, pm, l, health_list, health_of_1)
        current_mean_health = healthMean(N, health_list)
        print("Mean health after mutation:")
        print(current_mean_health)
        health_during_generations.append(current_mean_health)
        if np.abs(previous_mean_health - current_mean_health) < Properties.PRECISION and i > 1:
            counter = counter + 1
        else:
            counter = 0
        if counter >= 10:
            print("Mean health stays the same during 10 iterations -> algorithm was stopped")
            print("Iteration number: ", counter)
            build_line_graph(health_during_generations, N, l, x, y, type_of_selection, pm)
            save_to_file_the_end(1, 1, i, health_list, N, l, pm, type_of_selection)
            break
        if i == Properties.CONST_STOP_ALGORITHM:
            save_to_file_the_end(1, 1, i, health_list, N, l, pm, type_of_selection)
            build_histograms(pop, N, l, i, x, y, type_of_selection, pm, health_list)


run_genetic_algorithm_with_roulette(4, 100, 0, 1, mutation_probability)'''


def manage_patogenic_ch(i, pathogenic_muted_bool, border_patogenic_mation, pathogenic_muted_numb):
    if pathogenic_muted_bool:
        pathogenic_muted_numb += 1
        if pathogenic_muted_numb > border_patogenic_mation:
            return i


def run_genetic_algorithm_with_roulette(attempt, l, N, X, Y, pm, health_func, init_type, features=None):
    type_of_selection = 'Roulette'
    pop = generate_population(l, N, X, Y)
    counter = 0
    x = X * 100
    y = Y * 100
    mean_health_during_generations = []
    health_list = init_health_list(pop, health_func)

    border_patogenic_mutation = 0
    pathogenic_muted_numb = 0
    indexes_patogenic_muted = None

    if features:
        border_patogenic_mutation = num_of_pathogenic_genes(N, l, 0.0232)
        indexes_patogenic_muted = []

    for i in range(1, Properties.CONST_STOP_ALGORITHM + 1):

        previous_mean_health = healthMean(N, health_list)

        if i % Properties.CONST_FREQUENCY_PRINT_DIAGRAM == 0:
            build_histograms(pop, N, l, i, x, y, type_of_selection, pm, health_list, attempt, init_type)

        pop, health_list = selRoulette(pop, health_func, health_list, N)
        pop, health_list, pathogenic_muted_bool = mutation(pop, pm, l, health_list, health_func, features)

        if features:
            indexes_patogenic_muted.append(
                manage_patogenic_ch(i, pathogenic_muted_bool, border_patogenic_mutation, pathogenic_muted_numb))

        current_mean_health = healthMean(N, health_list)
        mean_health_during_generations.append(current_mean_health)

        if np.abs(previous_mean_health - current_mean_health) < Properties.PRECISION and i > 1:
            counter = counter + 1
        else:
            counter = 0
        if counter >= 10:
            build_histograms(pop, N, l, i, x, y, type_of_selection, pm, health_list, attempt, init_type)
            build_third_histogram(pop, N, l, i, x, y, type_of_selection, pm, attempt, init_type) # Distances to the wild type at the end of epoch
            build_line_graph(mean_health_during_generations, N, l, x, y, type_of_selection, pm, attempt, init_type)
            save_to_file_the_end(attempt, i, health_list, N, l, pm, type_of_selection, X, Y, -1,
                                 indexes_patogenic_muted)
            break
        if i == Properties.CONST_STOP_ALGORITHM:
            build_histograms(pop, N, l, i, x, y, type_of_selection, pm, health_list, attempt, init_type)
            build_third_histogram(pop, N, l, i, x, y, type_of_selection, pm, attempt, init_type)
            build_line_graph(mean_health_during_generations, N, l, x, y, type_of_selection, pm, attempt, init_type)
            save_to_file_the_end(attempt, i, health_list, N, l, pm, type_of_selection, X, Y, -1,
                                 indexes_patogenic_muted)

def run_genetic_algorithm_with_tournament(attempt, l, N, X, Y, pm, health_func, init_type, t, features=None):
    type_of_selection = 'Tournament'
    pop = generate_population(l, N, X, Y)
    counter = 0
    x = X * 100
    y = Y * 100
    mean_health_during_generations = []
    health_list = init_health_list(pop, health_func)

    border_patogenic_mutation = 0
    pathogenic_muted_numb = 0
    indexes_patogenic_muted = None

    if features:
        border_patogenic_mutation = num_of_pathogenic_genes(N, l, 0.0232)
        indexes_patogenic_muted = []

    for i in range(1, Properties.CONST_STOP_ALGORITHM + 1):
        previous_mean_health = healthMean(N, health_list)
        if i % Properties.CONST_FREQUENCY_PRINT_DIAGRAM == 0:
            build_histograms(pop, N, l, i, x, y, type_of_selection, pm, health_list, attempt, init_type)
        pop, health_list = tournament_selection(pop, health_func, health_list, t)
        pop, health_list, pathogenic_muted_bool = mutation(pop, pm, l, health_list, health_func, features)

        if features:
            indexes_patogenic_muted.append(
                manage_patogenic_ch(i, pathogenic_muted_bool, border_patogenic_mutation, pathogenic_muted_numb))

        current_mean_health = healthMean(N, health_list)
        mean_health_during_generations.append(current_mean_health)

        if np.abs(previous_mean_health - current_mean_health) < Properties.PRECISION and i > 1:
            counter = counter + 1
        else:
            counter = 0
        if counter >= 10:
            build_third_histogram(pop, N, l, i, x, y, type_of_selection, pm, attempt,
                                  init_type)  # Distances to the wild type at the end of epoch
            build_line_graph(mean_health_during_generations, N, l, x, y, type_of_selection, pm, attempt, init_type)
            save_to_file_the_end(attempt, i, health_list, N, l, pm, type_of_selection, X, Y, t,
                                 indexes_patogenic_muted)
            break
        if i == Properties.CONST_STOP_ALGORITHM:
            build_histograms(pop, N, l, i, x, y, type_of_selection, pm, health_list, attempt, init_type)
            build_third_histogram(pop, N, l, i, x, y, type_of_selection, pm, attempt, init_type)
            build_line_graph(mean_health_during_generations, N, l, x, y, type_of_selection, pm, attempt, init_type)
            save_to_file_the_end(attempt, i, health_list, N, l, pm, type_of_selection, X, Y, t,
                                 indexes_patogenic_muted)

def mutation_probabilities_for_roulette(l):
    px = 1 / (10 * l)
    #return [px + 0.2 * px, px - 0.2 * px, px / 2, px / 10, px / 100]
    return [0.000507, 0.001113, 0.000452, 0.000295] # constants from the file

def mutation_probabilities_for_tournament(t, l):
    if t == 2:
        return [0.001743, 0.000861, 0.00085, 0.000544]
    if t == 12:
        return [0.000574, 0.001749, 0.000881]
    if t == 4:
        px = 1 / (10 * l)
        return [px + 0.2 * px, px - 0.2 * px, px / 2, px / 10, px / 100]




# I do for selRoulette and tour = 12
def perform():
    N_1 = [100, 200, 800, 1000, 2000]
    N_2 = [100, 200, 800, 1000, 2000]
    l_N = [(10, N_1), (20, N_1), (80, N_1), (100, N_2), (200, N_2), (800, N_2), (1000, N_2), (2000, N_2), (8000, N_2)]

    # X and Y
    init_1 = [(0, 100)]
    init_2 = [(100, 0), (90, 10), (0, 100), (10, 90), (50, 50)]
    init_3 = [(100, 0), (90, 10), (50, 50)]

    pop_ratio = [init_1, init_2, init_3]
    health_funcs = [health_of_1, health_of_2, health_of_3]
    features = None

    for type_ind, list_ratios in enumerate(pop_ratio):
        for x, y in list_ratios:
            for l, list_N in l_N:
                arr_mutation_prob = mutation_probabilities_for_roulette(l) # ми ж хотіли використовувати контстанти з файлу
                for n in list_N:
                    if type_ind == 2:  # 3 type of init
                        features = list_features_of_ch(n, l)
                    for pm in arr_mutation_prob:
                        for attempt in range(3):
                            run_genetic_algorithm_with_roulette(attempt + 1, l, n, x / 100, y / 100, pm,
                                                                health_funcs[type_ind], type_ind, features)


#perform()
# Nika: I do for tournament t = 2 and t = 4

def perform_tournament():
    N_1 = [100, 200, 800, 1000, 2000]
    N_2 = [100, 200, 800, 1000, 2000]
    l_N = [(10, N_1), (20, N_1), (80, N_1), (100, N_2), (200, N_2), (800, N_2), (1000, N_2), (2000, N_2), (8000, N_2)]

    # X and Y
    init_1 = [(0, 100)]
    init_2 = [(100, 0), (90, 10), (0, 100), (10, 90), (50, 50)]
    init_3 = [(100, 0), (90, 10), (50, 50)]

    pop_ratio = [init_1, init_2, init_3]
    health_funcs = [health_of_1, health_of_2, health_of_3]
    features = None
    for t in [2, 4, 12]:
        for type_ind, list_ratios in enumerate(pop_ratio):
            for x, y in list_ratios:
                for l, list_N in l_N:
                    arr_mutation_prob = mutation_probabilities_for_tournament(t, l)
                    for n in list_N:
                        if type_ind == 2:  # 3 type of init
                            features = list_features_of_ch(n, l)
                        for pm in arr_mutation_prob:
                            for attempt in range(3):
                                run_genetic_algorithm_with_tournament(attempt + 1, l, n, x / 100, y / 100, pm,
                                                                      health_funcs[type_ind], type_ind, t, features)


perform_tournament()