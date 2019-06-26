import matplotlib.pyplot as plt
import numpy as np
import sys
import statistics
import csv
import math


def read_file(file):
    values = []
    csv_file = open(file)
    csv_reader = csv.reader(csv_file, delimiter=',')
    #  skipping first two lines line (containing metadata)
    meta = next(csv_reader)
    next(csv_reader)
    for row in csv_reader:
        values.append(row)
    print("reading done")
    return np.asarray(values, dtype=np.float32), np.asarray(meta, dtype=np.float32)


def convert_to_cluster(index):
    #  empty list for each cluster
    pos = -11 + index
    clusters = [[] for i in range(int(max(input[:, pos]) + 1))]
    for i in range(int(max(input[:, pos]) + 1)):
        clusters[i] = [s for s in input if s[pos] == i]

    #  converts clusters into a list of np-arrays and removes empty clusters
    delete_list = []
    for i in range(len(clusters)):
        if len(clusters[i]) == 0:
            delete_list.append(i)
        else:
            clusters[i] = np.asarray(clusters[i], dtype=np.float32)
    for i in delete_list:
        del clusters[i]

    #  print(clusters[0][:, pos])
    return clusters


######################################################################################################


#  for printing without dots
np.set_printoptions(threshold=sys.maxsize)
path = "./datasets_estimate/29.csv"
#  last entry is real class-label, 2. last entry is k-means class-label
input, meta_data = read_file(path)
print("meta data: " + str(meta_data))

number_of_features = int(meta_data[1])
perfect_number_clusters = int(meta_data[2])
meaningful_features = 5
rand_features = math.floor((number_of_features - meaningful_features) / 2)
const_features = math.ceil((number_of_features - meaningful_features) / 2)

functions = [np.amin, np.amax, np.ptp, np.var, np.std, np.mean, np.median, np.average]

print_list = []
for k_value in range(10):
    clusters = convert_to_cluster(k_value)

    result_clusters = [np.zeros(shape=(len(functions), number_of_features)) for i in range(len(clusters))]
    mean_result_clusters = [np.zeros(shape=(len(functions), 3)) for i in range(len(clusters))]
    result_dataset = np.zeros(shape=(len(functions), number_of_features))

    ind = 0
    for f in functions:
        for i in range(number_of_features):
            result_dataset[ind, i] = f(input[:, i])
        ind = ind + 1

    # print(f.__name__ + ": " + str(f(input[:, i])))
    for j in range(len(clusters)):
        ind = 0
        for f in functions:
            for i in range(number_of_features):
                result_clusters[j][ind, i] = f(clusters[j][:, i])
            ind = ind + 1
    #  print(result_dataset)
    #  print(result_clusters[3])

    for j in range(len(clusters)):
        ind = 0
        for f in functions:
            mean_result_clusters[j][ind, 0] = np.mean(result_clusters[j][ind, : meaningful_features])
            mean_result_clusters[j][ind, 1] = np.mean(
                result_clusters[j][ind, meaningful_features: meaningful_features + rand_features])
            mean_result_clusters[j][ind, 2] = np.mean(result_clusters[j][ind, meaningful_features + rand_features:])
            ind = ind + 1
    print_list.append(mean_result_clusters[0][3, 0])

figure = plt.figure()
ax1 = figure.add_subplot(111)
title = "Cluster 0"
ax1.set_title(title)
ax1.set_xlabel('Metrics')
ax1.set_ylabel('Metric-Values')
#  plots x-axis descending sorted


x_axis = [f.__name__ for f in functions]
y_axis = mean_result_clusters[0]
ax1.plot(x_axis, y_axis[:, 0], label="Meaningful Feature")
ax1.plot(x_axis, y_axis[:, 1], label="Random Feature")
ax1.plot(x_axis, y_axis[:, 2], label="Constant Feature")


# ax1.plot(x_axis, y_axis[2], label='Cluster 2')
ax1.legend(loc=1, fontsize='small')
# show only every 5th entry of x-axis
# r = np.arange(0, len(sorted_numbers), step=5)
# plt.xticks(r, [to_sort[i] for i in r], rotation=45)
# log-scale
# ax1.set_yscale('log')
# axis-range
# ax1.set_ylim([0.0001, 25])
png_name = "./plot/meanMetrikVerlauf.png"
plt.savefig(png_name)
plt.close(figure)


figure2 = plt.figure()
ax2 = figure2.add_subplot(111)
title = "Cluster 0"
ax2.set_title(title)
ax2.set_xlabel('Metrics')
ax2.set_ylabel('Metric-Values')
x_axis = [f.__name__ for f in functions]
y_axis = result_clusters[0]
for r in range(len(y_axis[0])):
    ax2.plot(x_axis, y_axis[:, r], label="Feature " + str(r))
# ax1.plot(x_axis, y_axis[2], label='Cluster 2')
ax2.legend(loc=1, fontsize='small')
png_name = "./plot/allMetrikVerlauf.png"
plt.savefig(png_name)
plt.close(figure2)


figure3 = plt.figure()
ax3 = figure3.add_subplot(111)
title = "Cluster 0"
ax3.set_title(title)
ax3.set_xlabel('Metrics')
ax3.set_ylabel('Metric-Values')
x_axis = [i + 1 for i in range(10)]
y_axis = print_list
ax3.plot(x_axis, y_axis, label="Variance")
ax3.legend(loc=1, fontsize='small')
png_name = "./plot/verlaufVarianzDifferentK.png"
plt.savefig(png_name)
plt.close(figure3)
