import numpy as np  # linear algebra
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
import matplotlib.pyplot as plt

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

'''Comparing two binary strings of equal length, Hamming distance is the number of bit positions in which the two bits are different.
In order to calculate the Hamming distance between two strings, we perform their XOR operation and then count the total number of 1s in the resultant string.
'''
def calculate_chromosome_fitness(chromosome, fitness_func="hamming distance"):
    if fitness_func=="hamming distance":
        pass
    else:
        pass

def mutation(N, l):
    px = 1/(10*l)
    pm = px
    #pm = px + 0.2*px
    #pm = px - 0.2*px
    #pm = px/2
    #pm = px/10
    #pm = px/100

def tournament_selection(t):
    if t == 2:
        pass
    elif t == 4:
        pass
    elif t == 12:
        pass

# testing
pop = generate_population(10, 100, 0.5, 0.5)
print(pop[:-20])

