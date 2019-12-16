import numpy as np

def generate_population(l, N, X, Y, coding="binary"):
    if coding == "binary":
        x = int(X * N)
        y = int(Y * N)
        population_with_zeros = np.random.randint(0, 1, (x, l))
        uniform_population = np.random.randint(0, 2, (y, l))
        general_population = np.concatenate((population_with_zeros, uniform_population), axis=0)
    return general_population


def health_of_1(ch):
    return np.count_nonzero(ch == 0)


pop = generate_population(10, 100, 0, 1)
health_list = [health_of_1(ch) for ch in pop]

health_list, pop = (list(t) for t in zip(*sorted(zip(health_list, pop))))
#print(pop, health_list)


