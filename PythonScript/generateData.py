import numpy as np
import math
from sklearn.cluster import k_means
from scipy.stats import cauchy


def generate_point(distribution, noise, cluster_number, number_of_clusters):
    max_value = number_of_clusters * 10
    overlap = (noise * np.random.uniform(low=0, high=1))
    min = (10*cluster_number) - (overlap*10)
    max = 10 + (overlap*10) + (10*cluster_number)
    if distribution == 'g':
        point = [np.random.normal(loc=(max+min)/2, scale=(max-min)/3),
                 np.random.normal(loc=(max+min)/2),
                 np.random.normal(loc=(max+min)/2, scale=(max-min)/3),
                 np.random.randint(max_value),
                 np.random.randint(4) * (max_value/4),
                 np.random.randint(2) * max_value,
                 12345,
                 cluster_number * 10]
    elif distribution == 'c':
        point = [cauchy.rvs(loc=(max+min)/2, scale=(max-min)/100),
                 cauchy.rvs(loc=(max + min) / 2, scale=(max - min) / 100),
                 cauchy.rvs(loc=(max + min) / 2, scale=(max - min) / 100),
                 np.random.randint(max_value),
                 np.random.randint(4) * (max_value / 4),
                 np.random.randint(2) * max_value,
                 12345,
                 cluster_number * 10]
    else:
        point = [np.random.uniform(low=(10*cluster_number), high=(10*cluster_number)+10),
                 np.random.uniform(low=min, high=max),
                 np.random.uniform(low=min, high=max),
                 np.random.randint(max_value),
                 np.random.randint(4) * (max_value / 4),
                 np.random.randint(2) * max_value,
                 12345,
                 cluster_number * 10]
    return point


def generate_data(distribution, noise, number_of_clusters, points_per_cluster):
    out_array = [generate_point(distribution, noise, math.floor(i/points_per_cluster), number_of_clusters)
                 for i in range(points_per_cluster * number_of_clusters)]
    return out_array


data = generate_data(distribution='c', noise=0, number_of_clusters=3, points_per_cluster=10)
for x in data:
    print(x)
kmeans = k_means(init='k-means++', X=data, n_clusters=3)
print(kmeans[1])
# appends k_means cluster label to each point
for i in range(len(data)):
    data[i].append(kmeans[1][i])
