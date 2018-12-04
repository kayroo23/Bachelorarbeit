import numpy as np
import math
from sklearn.cluster import k_means


def generate_point(distribution, noise, cluster_number, number_of_clusters):
    max_value = number_of_clusters * 10
    if distribution == 'g':
        point = [np.random.normal(), np.random.normal(), np.random.normal(), np.random.randint(max_value),
                 np.random.randint(4) * (max_value/4), np.random.randint(2) * max_value, 12345, cluster_number * 10]
    elif distribution == 'c':
        point = [np.random.standard_cauchy(), np.random.normal()]
    else:
        point = [np.random.uniform(), np.random.normal()]

    return point


def generate_data(distribution, noise, number_of_clusters, points_per_cluster):
    out_array = [generate_point(distribution, noise, math.floor(i/points_per_cluster), number_of_clusters) for i in range(points_per_cluster * number_of_clusters)]
    return out_array


data = generate_data('g', 0, 2, 2)
print(data)
kmeans = k_means(X=data, n_clusters=2)
print(kmeans)
