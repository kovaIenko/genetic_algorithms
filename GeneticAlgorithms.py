import random
import Properties
import csv
import numpy as np  # linear algebra
from collections import Counter
from decimal import *
import os
import matplotlib.pyplot as plt
from math import floor

getcontext().prec = 6

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
def health_of_3(ch, type_=None):
    l = len(ch)
    if type_ == Properties.CONST_TYPE_NEUTRAL:
        return l
    elif type_ == Properties.CONST_TYPE_LETHAL:
        return 0.1
    else:
        return l - l * Properties.CONST_K / 100

def health_of_2(ch):
    l = len(ch)
    return l

def sort_sort(pop, list_of_health):
    if list_of_health:
        list_of_health, pop = zip(*sorted(zip(list_of_health, pop), key=lambda x: x[0]))
        return list(pop), list(list_of_health)
    else:
        return pop, list_of_health

# RWS
def selRoulette(pop, list_of_health, N):
    chosen = []
    pop, list_of_health = sort_sort(pop, list_of_health)
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


def mutation(pop, percentage, l, list_of_health=None, features=None):
    count_ = len(pop) * l
    numbOfMut = int(round(count_ * percentage))
    pathogenic_muted_counter = 0
    for i in range(numbOfMut):
        index_of_ch, index_of_gen = generate_ch_and_gen(count_, l)
        if not features and list_of_health:
            currCh = pop[index_of_ch]
            currentCh = turnOverGen(currCh, index_of_gen)
            pop[index_of_ch] = currentCh
            list_of_health[index_of_ch] = health_of_1(currentCh)
        elif not list_of_health:
            currCh = pop[index_of_ch]
            currentCh = turnOverGen(currCh, index_of_gen)
            pop[index_of_ch] = currentCh
        else:
            current_map = features[index_of_ch]
            type_ = getType(current_map, index_of_gen)
            pathogenic_muted_counter = increment_neutral_counter(pathogenic_muted_counter, type_)
            list_of_health[index_of_ch] = health_of_3(type_, l)
    return pop, list_of_health, pathogenic_muted_counter


def increment_neutral_counter(counter, type):
    if type == Properties.CONST_TYPE_NEUTRAL:
        counter += 1
    return counter


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
    for type_, list_of_genes in current_map.items():
        if ind_of_gen in list_of_genes:
            return type_
    return None


# list of maps for each chromosome
def list_features_of_ch(N, l):
    maps = []
    for ch in range(N):
        maps.append(mutation_genes_distribution(l))
    return maps


# Returns a map for a chromosome, where key is a mutation type and value is a list of genes positions
def mutation_genes_distribution(l, first_neutral_percent=0.135, any_neutral_percent=0.245,
                                pathogenic_percent=Properties.CONST_PATHOGENIC_PERCENT):
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
def calc_all_distances(pop, N_, l_):
    frequency = {new_list: 0 for new_list in range(l_ + 1)}
    for i in range(0, N_):
        for j in range(i + 1, N_):
            dis = hamming(pop[i], pop[j])
            frequency[dis] = frequency.get(dis) + 1
    return frequency


def healthMean(N, l, list_health):
    if not list_health:
        return l
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
                               'HammingHist;InitType={0};Attempt={1};Selection={2},N={3};l={4}/'.format(init_type,
                                                                                                        attempt,
                                                                                                        selection_type,
                                                                                     N, l))
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    distances = calc_all_distances(pop, N, l)
    plt.bar(list(distances.keys()), distances.values(), color='g', width=0.9)
    plt.xticks(list(distances.keys()))
    plt.savefig(results_dir + "/X={0};Y={1};iter={2};pm={3}.png".format(x, y, iter_num, pm))
    plt.clf()


# Save a histogram to png file with all the parameters specified
# (a histogram with frequencies of Hamming distances to the optimal chromosome from each chromosome in the population)

def build_second_histogram(pop, list_health, N, l, iter_num, x, y, selection_type, pm, attempt, init_type):
    script_dir = os.path.dirname(__file__) + '/Plots'
    results_dir = os.path.join(script_dir,
                               'OptimalHist;InitType={0};Attempt={1};Selection={2},N={3};l={4}/'.format(init_type,
                                                                                                        attempt,
                                                                                                        selection_type,
                                                                                                        N, l))
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    distances = calc_hamming_to_ideal(pop, list_health, l, init_type)
    plt.bar(list(distances.keys()), distances.values(), color='red', width=0.9)
    plt.xticks(list(distances.keys()))
    plt.savefig(results_dir + "/X={0};Y={1};iter={2};pm={3}.png".format(x, y, iter_num, pm))
    plt.clf()


# Save a histogram to png file with all the parameters specified
# (a histogram with frequencies of Hamming distances to the optimal chromosome from each chromosome in the population)

def build_third_histogram(pop, N, l, iter_num, x, y, selection_type, pm, attempt, init_type):
    script_dir = os.path.dirname(__file__) + '/Plots'
    results_dir = os.path.join(script_dir,
                               'WildTypeHist;InitType={0};Attempt={1};Selection={2},N={3};l={4}/'.format(init_type,
                                                                                                         attempt,
                                                                                                         selection_type,
                                                                                                         N, l))
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
                               'LineGraph;InitType={0};Attempt={1};Selection={2},N={3};l={4};/'.format(init_type,
                                                                                                       attempt,
                                                                                                       selection_type,
                                                                                                       N, l))
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


def init_health_list(pop, l, N, method=1):
    list_of_health = []
    if method == 1:
        for ch in pop:
            list_of_health.append(health_of_1(ch))
        return list_of_health
    elif method == 2:
        return None
    else:
        return [l] * N


def encode_key(numb):
    return int(numb * 10000)


def calc_hamming_to_ideal(pop, list_health, l, method=1):
    if method == 1:
        frequency = {new_list: 0 for new_list in range(l + 1)}
        for ch_health in list_health:
            dis = l - ch_health
            frequency[dis] = frequency.get(dis) + 1
    elif method == 2:
        frequency = {new_list: 0 for new_list in range(l + 1)}
        for ch in pop:
            dis = l - health_of_1(ch)
            frequency[dis] = frequency.get(dis) + 1
    else:
        frequency = {0: 0, encode_key(l - 0.1): 0, encode_key(l - (Decimal(l - Properties.CONST_K * l / 100))): 0}
        for ch_health in list_health:
            dis = l - ch_health
            value_ = frequency.get(encode_key(dis))
            frequency.update({encode_key(dis): 1 if not value_ else value_ + 1})
        frequency = dict((float(Decimal(key) / Decimal(10000)), value) for (key, value) in frequency.items())
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


def bestHealth(list_of_health, l):
    if not list_of_health:
        return l
    return Decimal(max(list_of_health))


def deviation_meanHealth_and_optimum(list_of_health, N, l):
    mean = healthMean(N, l, list_of_health)
    return Decimal(l) - mean


def deviation_bestHealth_and_optimum(list_of_health, l):
    return l - bestHealth(list_of_health, l)


def percent_polym_genes(pop, list_of_health, l, N):
    if not list_of_health:
        list_of_health = init_health_list(pop, l, N, 1)
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


def save_to_file_the_end(pop, attempt, iteration, list_of_health, N, l, pm, type_of_selection, X, Y, init_type,
                         params_tour=-1, arr_neutral_iter=None):
    # Another variant with dynamically defining csv columns

    csv_columns = ['attempt', 'method', 'iteration', 'N', 'l', 'X', 'Y', 'pm', 'selection', 'health_mean', 'health_best', 'deviation_meanHealth_and_optimum', 'deviation_bestHealth_and_optimum', 'percent_polym_genes', 'arr_neutral_iter']
    row_map = {"attempt": attempt, "method": init_type, "N": N, "l": l, "X": X, "Y": Y, "iteration": iteration, "pm": pm}
    if type_of_selection == "Tournament":
        type_of_selection = type_of_selection + ' ' + str(params_tour)
    row_map["selection"] = type_of_selection
    if arr_neutral_iter:
        row_map["arr_neutral_iter"] = arr_neutral_iter
    row_map["health_mean"] = healthMean(N, l, list_of_health)
    row_map["health_best"] = bestHealth(list_of_health, l)
    row_map["deviation_meanHealth_and_optimum"] = deviation_meanHealth_and_optimum(list_of_health, N, l)
    row_map["deviation_bestHealth_and_optimum"] = deviation_bestHealth_and_optimum(list_of_health, l)
    row_map["percent_polym_genes"] = percent_polym_genes(pop, list_of_health, l, N)
    csv_file = "data.csv"
    try:
        with open(csv_file, 'a') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
            if os.stat(csv_file).st_size == 0:
                writer.writeheader()
            writer.writerow(row_map)
    except IOError:
        print("I/O error")

    print(row_map)



def should_be_stopped(i, sum_mean_health, current_mean_health):
    if i % Properties.CONST_NUMB_GETTING_INFO == 0 and i != 0:  # each 11n iterations
        last_10_mean_health = sum_mean_health / Properties.CONST_NUMB_GETTING_INFO
        if Decimal(current_mean_health) - Decimal(last_10_mean_health) < Properties.PRECISION:
            print("Algorithm was stopped: there is mistake ")
            return True
    return False


# TO CHANGE TO NEUTRAL
# how many mutations we can do through the all iterations
def num_of_neutral_genes(N, l, percentage):
    return N * l * percentage


# At each h iterations we build only histograms for pairwise Hamming distances & Hamming distances to the optimal
# chromosome

def build_histograms(pop, N, l, i, x, y, type_of_selection, pm, health_list, attempt, init_type):
    build_first_histogram(pop, N, l, i, x, y, type_of_selection, pm, attempt, init_type)
    build_second_histogram(pop, health_list, N, l, i, x, y, type_of_selection, pm, attempt, init_type)


def manage_pathogenic_ch(i, arr_of_indexes, neutral_muted_counter, limit_neutral_mutation):
    if neutral_muted_counter > limit_neutral_mutation:
        arr_of_indexes.append(i)
        neutral_muted_counter = 0
    return arr_of_indexes, neutral_muted_counter


def run_genetic_algorithm_with_roulette(attempt, l, N, x, y, pm, init_type, features=None):
    type_of_selection = 'Roulette'
    pop = generate_population(l, N, x/100, y/100)
    counter = 0
    mean_health_during_generations = []
    health_list = init_health_list(pop, l, N, init_type)

    border_neutral_mutation = 0
    muted_counter = 0
    arr_neutral_indexes = None

    if features:
        border_neutral_mutation = num_of_neutral_genes(N, l, Properties.CONST_NEUTRAL_PERCENT)
        arr_neutral_indexes = []

    for i in range(1, Properties.CONST_STOP_ALGORITHM + 1):

        previous_mean_health = healthMean(N, l, health_list) if init_type != 2 else l

        if i % Properties.CONST_FREQUENCY_PRINT_DIAGRAM == 0:
            build_histograms(pop, N, l, i, x, y, type_of_selection, pm, health_list, attempt, init_type)

        pop, health_list = selRoulette(pop, health_list, N)
        pop, health_list, neutral_muted_counter = mutation(pop, pm, l, health_list, features)

        if features:
            muted_counter += neutral_muted_counter
            arr_neutral_indexes, muted_counter = manage_pathogenic_ch(i, arr_neutral_indexes,
                                                                              muted_counter,
                                                                              border_neutral_mutation)
        current_mean_health = healthMean(N, l, health_list)
        mean_health_during_generations.append(current_mean_health)
        if init_type != 2:  # we don't need to stop algorithm because of similar mean health in case of the 2nd init type
            if np.abs(previous_mean_health - current_mean_health) < Properties.PRECISION and i > 1:
                counter = counter + 1
            else:
                counter = 0
        # If mean health between populations doesn't differ much
        if counter >= 10:
            build_histograms(pop, N, l, i, x, y, type_of_selection, pm, health_list, attempt, init_type)
            build_third_histogram(pop, N, l, i, x, y, type_of_selection, pm, attempt,
                                  init_type)  # Distances to the wild type at the end of epoch
            build_line_graph(mean_health_during_generations, N, l, x, y, type_of_selection, pm, attempt, init_type)
            save_to_file_the_end(pop, attempt, i, health_list, N, l, pm, type_of_selection, x, y, init_type, -1,
                                 arr_neutral_indexes)
            break
        # If final iteration
        if i == Properties.CONST_STOP_ALGORITHM:
            build_histograms(pop, N, l, i, x, y, type_of_selection, pm, health_list, attempt, init_type)
            build_third_histogram(pop, N, l, i, x, y, type_of_selection, pm, attempt, init_type)
            build_line_graph(mean_health_during_generations, N, l, x, y, type_of_selection, pm, attempt, init_type)
            save_to_file_the_end(pop, attempt, i, health_list, N, l, pm, type_of_selection, x, y, init_type,  -1,
                                 arr_neutral_indexes)

def mutation_probabilities_for_roulette(l):
    px = 1 / (10 * l)
    # return [px + 0.2 * px, px - 0.2 * px, px / 2, px / 10, px / 100]
    return [0.000507, 0.001113, 0.000452, 0.000295]  # constants from the file


def mutation_probabilities_for_tournament(t, l):
    if t == 2:
        return [0.001743, 0.000861, 0.00085, 0.000544]
    if t == 12:
        return [0.000574, 0.001749, 0.000881]
    if t == 4:
        px = 1 / (10 * l)
        return [px + 0.2 * px, px - 0.2 * px, px / 2, px / 10, px / 100]


def get_init_data():
    N_1 = [100, 200]  # , 800, 1000, 2000]
    N_2 = [100, 200]  # , 800, 1000, 2000]
    l_N = [(10, N_1), (20, N_1), (80, N_1), (100, N_2),
           (200, N_2)]  # , (800, N_2), (1000, N_2), (2000, N_2), (8000, N_2)]

    # X and Y
    init_1 = [(0, 100)]
    init_2 = [(100, 0), (90, 10), (0, 100), (10, 90), (50, 50)]
    init_3 = [(100, 0), (90, 10), (50, 50)]

    pop_ratio = [init_1, init_2, init_3]
    return l_N, pop_ratio


# I do for selRoulette and tour = 12
def perform_roulette():
    l_N, pop_ratio = get_init_data()
    features = None

    for type_ind, list_ratios in enumerate(pop_ratio):
        for x, y in list_ratios:
            for l, list_N in l_N:
                arr_mutation_prob = mutation_probabilities_for_roulette(
                    l)  # ми ж хотіли використовувати контстанти з файлу
                for n in list_N:
                    if type_ind == 2:  # 3 type of init
                        features = list_features_of_ch(n, l)
                    for pm in arr_mutation_prob:
                        for attempt in range(3):
                          run_genetic_algorithm_with_roulette(attempt + 1, l, n, x, y , pm, type_ind+1, features)

perform_roulette()
