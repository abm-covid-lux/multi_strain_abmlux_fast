import csv
import os
import numpy as np
import math
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

def plot_R_differences(inf_multipliers, trans_multipliers, R_differences, show_fig, save_fig):
    """Plots R_differences as a heat map"""

    fig, ax = plt.subplots()
    fig.set_size_inches(16, 16)
    cmap = plt.get_cmap('PuOr') # RdYlGn
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

    R_differences = np.rot90(R_differences)

    bound = max(abs(np.min(R_differences)), abs(np.max(R_differences)))
    bound = 150000
    bounds = np.arange(-bound,bound, math.ceil((2*bound + 1)/256))
    idx=np.searchsorted(bounds,0)
    bounds=np.insert(bounds,idx,0)
    norm = BoundaryNorm(bounds, cmap.N)

    trans_multipliers.reverse()

    ax.set_xticks(np.arange(len(inf_multipliers)))
    ax.set_yticks(np.arange(len(trans_multipliers)))
    ax.set_xticklabels(inf_multipliers, fontsize=18)
    ax.set_yticklabels(trans_multipliers, fontsize=18)

    plt.setp(ax.get_xticklabels(), rotation=90, ha="center", rotation_mode="default")
    plt.xlabel('Infectious period ratio', fontsize=18)
    plt.ylabel('Transmission probability ratio', fontsize=18)
    plt.title('Reduction in Total Cumulative Cases', fontsize=18)
    plt.imshow(R_differences,interpolation='none',norm=norm,cmap=cmap)
    cbar = plt.colorbar(ticks=np.linspace(-bound, bound, 11), fraction=0.04525)
    cbar.ax.tick_params(labelsize=18)
    if show_fig:
        plt.show()
    if save_fig:
        fig.savefig('differences.png', bbox_inches='tight')

def get_cumulative_count(data):
    """Extract cumulative count data for two strains"""

    cumulative_count = int(data[-1][0]) + int(data[-1][1])

    return cumulative_count

cumulative_counts = defaultdict(lambda: defaultdict(list))

# output_directories = ['output_prng_1', 'output_prng_2', 'output_prng_3', 'output_prng_4',\
#                       'output_prng_5', 'output_prng_6', 'output_prng_7', 'output_prng_8',\
#                       'output_prng_9']

output_directories = ['output/output_prng_1']

for output_directory in output_directories:
    for filename in os.listdir(output_directory):
        f = os.path.join(output_directory, filename)
        if os.path.isfile(f):
            if filename[0:27] == "resident_cumulative_counts_":
                with open(f, newline='') as csvfile:
                    cumulative_counts[filename[38: -11]][filename[27:38]].append(\
                        get_cumulative_count(list(csv.reader(csvfile))))

differences = defaultdict(dict)
for param_key in cumulative_counts:
    num_scenario_0 = len(cumulative_counts[param_key]["scenario_0_"])
    num_scenario_1 = len(cumulative_counts[param_key]["scenario_1_"])
    difference = (sum(cumulative_counts[param_key]["scenario_0_"])/num_scenario_0) -\
                 (sum(cumulative_counts[param_key]["scenario_1_"])/num_scenario_1)
    trans_multiplier = float(param_key[14:17])
    inf_multiplier = float(param_key[4:7])
    differences[inf_multiplier][trans_multiplier] = float(difference)

inf_multipliers = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4,
                   2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8]
trans_multipliers = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4,
                     2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8]

R_differences = []
for inf_multiplier in inf_multipliers:
    R_differences_inf_multiplier = []
    for trans_multiplier in trans_multipliers:
        if inf_multiplier in differences and trans_multiplier in differences[inf_multiplier]:
            R_differences_inf_multiplier.append(differences[inf_multiplier][trans_multiplier])
        else:
            R_differences_inf_multiplier.append(0)
    R_differences.append(R_differences_inf_multiplier)

plot_R_differences(inf_multipliers, trans_multipliers, np.array(R_differences), False, True)
