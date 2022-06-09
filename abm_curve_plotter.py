import csv
import os
import numpy as np
from collections import defaultdict
from tqdm import tqdm
import matplotlib.pyplot as plt

def get_strain_counts(data):
    """Extract strain count data for two strains"""

    strain_counts_0 = []
    strain_counts_1 = []

    for row in data[1:]:
        strain_counts_0.append(int(row[0]))
        strain_counts_1.append(int(row[1]))

    strain_counts = [strain_counts_0, strain_counts_1]

    return np.array(strain_counts)

# Figure 7 (left)
# baseline = 'Baseline'
# intervention = 'Lockdown'
# output_baseline = 'resident_strain_counts_scenario_0'
# output_intervention = 'resident_strain_counts_scenario_1'
# baseline_directory = 'output/example_1'
# intervention_directory = 'output/example_1'
# strain_counts = defaultdict(dict)
# for filename in os.listdir(baseline_directory):
#     f = os.path.join(baseline_directory, filename)
#     if os.path.isfile(f):
#         if filename[0:len(output_baseline)] == output_baseline:
#             with open(f, newline='') as csvfile:
#                 strain_counts[baseline][f] = get_strain_counts(list(csv.reader(csvfile)))
# for filename in os.listdir(intervention_directory):
#     f = os.path.join(intervention_directory, filename)
#     if os.path.isfile(f):
#         if filename[0:len(output_intervention)] == output_intervention:
#             with open(f, newline='') as csvfile:
#                 strain_counts[intervention][f] = get_strain_counts(list(csv.reader(csvfile)))
# T = max([len(strain_counts[baseline][filename][0]) for filename in strain_counts[baseline]])
# strain_counts_averages = {baseline: np.average([strain_counts[baseline][filename] for filename in strain_counts[baseline]], axis=0),
#                           intervention: np.average([strain_counts[intervention][filename] for filename in strain_counts[intervention]], axis=0)}
# labels = list(range(T))
# strain_1 = 'Recovered'
# strain_2 = 'Recovered'
# plt.figure(figsize=(8, 6))
# font = {'size' : 12}
# plt.rc('font', **font)
# x_lim = 300*144
# plt.xlim(0, x_lim)
# plt.ylim(0, 40000)
# for filename in strain_counts[baseline]:
#     plt.plot(list(range(T)), np.add(strain_counts[baseline][filename][0], strain_counts[baseline][filename][1]), "black", linewidth=1, alpha=0.1)
# plt.plot(list(range(T)), np.add(strain_counts_averages[baseline][0], strain_counts_averages[baseline][1]), "black", linewidth=1, alpha=1.0, label='Infected: ' + baseline)
# for filename in strain_counts[intervention]:
#     plt.plot(list(range(T)), np.add(strain_counts[intervention][filename][0], strain_counts[intervention][filename][1]), "black", linestyle='dotted', linewidth=1, alpha=0.1)
# plt.plot(list(range(T)), np.add(strain_counts_averages[intervention][0], strain_counts_averages[intervention][1]), "black", linestyle='dotted', linewidth=1, alpha=1.0, label='Infected: ' + intervention)
# plt.xlabel('Day')
# plt.xticks(ticks=[a*28*144 for a in range((x_lim // (28*144)) + 1)],
#            labels=[a*28 for a in range((x_lim // (28*144)) + 1)])
# plt.legend(loc='upper right')
# plt.grid(False)
# plt.savefig('fig7left.png', bbox_inches='tight')

# Figure 7 (right)
# baseline = 'Baseline'
# intervention = 'Lockdown'
# output_baseline = 'resident_cumulative_counts_scenario_0'
# output_intervention = 'resident_cumulative_counts_scenario_1'
# baseline_directory = 'output/example_1'
# intervention_directory = 'output/example_1'
# strain_counts = defaultdict(dict)
# for filename in os.listdir(baseline_directory):
#     f = os.path.join(baseline_directory, filename)
#     if os.path.isfile(f):
#         if filename[0:len(output_baseline)] == output_baseline:
#             with open(f, newline='') as csvfile:
#                 strain_counts[baseline][f] = get_strain_counts(list(csv.reader(csvfile)))
# for filename in os.listdir(intervention_directory):
#     f = os.path.join(intervention_directory, filename)
#     if os.path.isfile(f):
#         if filename[0:len(output_intervention)] == output_intervention:
#             with open(f, newline='') as csvfile:
#                 strain_counts[intervention][f] = get_strain_counts(list(csv.reader(csvfile)))
# T = max([len(strain_counts[baseline][filename][0]) for filename in strain_counts[baseline]])
# strain_counts_averages = {baseline: np.average([strain_counts[baseline][filename] for filename in strain_counts[baseline]], axis=0),
#                           intervention: np.average([strain_counts[intervention][filename] for filename in strain_counts[intervention]], axis=0)}
# labels = list(range(T))
# strain_1 = 'Recovered'
# strain_2 = 'Recovered'
# plt.figure(figsize=(8, 6))
# font = {'size' : 12}
# plt.rc('font', **font)
# x_lim = 300*144
# plt.xlim(0, x_lim)
# plt.ylim(0, 200000)
# for filename in strain_counts[baseline]:
#     plt.plot(list(range(T)), np.add(strain_counts[baseline][filename][0], strain_counts[baseline][filename][1]), "black", linewidth=1, alpha=0.1)
# plt.plot(list(range(T)), np.add(strain_counts_averages[baseline][0], strain_counts_averages[baseline][1]), "black", linewidth=1, alpha=1.0, label='Recovered: ' + baseline)
# for filename in strain_counts[intervention]:
#     plt.plot(list(range(T)), np.add(strain_counts[intervention][filename][0], strain_counts[intervention][filename][1]), "black", linestyle='dotted', linewidth=1, alpha=0.1)
# plt.plot(list(range(T)), np.add(strain_counts_averages[intervention][0], strain_counts_averages[intervention][1]), "black", linestyle='dotted', linewidth=1, alpha=1.0, label='Recovered: ' + intervention)
# plt.xlabel('Day')
# plt.xticks(ticks=[a*28*144 for a in range((x_lim // (28*144)) + 1)],
#            labels=[a*28 for a in range((x_lim // (28*144)) + 1)])
# plt.legend(loc='upper right')
# plt.grid(False)
# plt.savefig('fig7right.png', bbox_inches='tight')

# Figure 9 (left)
# baseline = 'Baseline'
# intervention = 'Lockdown'
# output_baseline = 'resident_strain_counts_scenario_0'
# output_intervention = 'resident_strain_counts_scenario_1'
# baseline_directory = 'output/example_2'
# intervention_directory = 'output/example_2'
# strain_counts = defaultdict(dict)
# for filename in os.listdir(baseline_directory):
#     f = os.path.join(baseline_directory, filename)
#     if os.path.isfile(f):
#         if filename[0:len(output_baseline)] == output_baseline:
#             with open(f, newline='') as csvfile:
#                 strain_counts[baseline][f] = get_strain_counts(list(csv.reader(csvfile)))
# for filename in os.listdir(intervention_directory):
#     f = os.path.join(intervention_directory, filename)
#     if os.path.isfile(f):
#         if filename[0:len(output_intervention)] == output_intervention:
#             with open(f, newline='') as csvfile:
#                 strain_counts[intervention][f] = get_strain_counts(list(csv.reader(csvfile)))
# T = max([len(strain_counts[baseline][filename][0]) for filename in strain_counts[baseline]])
# strain_counts_averages = {baseline: np.average([strain_counts[baseline][filename] for filename in strain_counts[baseline]], axis=0),
#                           intervention: np.average([strain_counts[intervention][filename] for filename in strain_counts[intervention]], axis=0)}
# labels = list(range(T))
# strain_1 = 'Strain 1'
# strain_2 = 'Strain 2'
# plt.figure(figsize=(8, 6))
# font = {'size' : 12}
# plt.rc('font', **font)
# x_lim = 300*144
# plt.xlim(0, x_lim)
# plt.ylim(0, 80000)
# for filename in strain_counts[baseline]:
#     plt.plot(list(range(T)), strain_counts[baseline][filename][0], "red", linewidth=1, alpha=0.1)
#     plt.plot(list(range(T)), strain_counts[baseline][filename][1], "blue", linewidth=1, alpha=0.1)
# plt.plot(list(range(T)), strain_counts_averages[baseline][0], "red", linewidth=1, alpha=1.0, label='Infected: ' + baseline + ': Strain 1')
# plt.plot(list(range(T)), strain_counts_averages[baseline][1], "blue", linewidth=1, alpha=1.0, label='Infected: ' + baseline + ': Strain 2')
# for filename in strain_counts[intervention]:
#     plt.plot(list(range(T)), strain_counts[intervention][filename][0], "red", linestyle='dotted', linewidth=1, alpha=0.1)
#     plt.plot(list(range(T)), strain_counts[intervention][filename][1], "blue", linestyle='dotted', linewidth=1, alpha=0.1)
# plt.plot(list(range(T)), strain_counts_averages[intervention][0], "red", linestyle='dotted', linewidth=1, alpha=1.0, label='Infected: ' + intervention + ': Strain 1')
# plt.plot(list(range(T)), strain_counts_averages[intervention][1], "blue", linestyle='dotted', linewidth=1, alpha=1.0, label='Infected: ' + intervention + ': Strain 2')
# plt.xlabel('Day')
# plt.xticks(ticks=[a*28*144 for a in range((x_lim // (28*144)) + 1)],
#            labels=[a*28 for a in range((x_lim // (28*144)) + 1)])
# plt.legend(loc='upper right')
# plt.grid(False)
# plt.savefig('fig9left.png', bbox_inches='tight')

# Figure 9 (right)
# baseline = 'Baseline'
# intervention = 'Lockdown'
# output_baseline = 'resident_cumulative_counts_scenario_0'
# output_intervention = 'resident_cumulative_counts_scenario_1'
# baseline_directory = 'output/example_2'
# intervention_directory = 'output/example_2'
# strain_counts = defaultdict(dict)
# for filename in os.listdir(baseline_directory):
#     f = os.path.join(baseline_directory, filename)
#     if os.path.isfile(f):
#         if filename[0:len(output_baseline)] == output_baseline:
#             with open(f, newline='') as csvfile:
#                 strain_counts[baseline][f] = get_strain_counts(list(csv.reader(csvfile)))
# for filename in os.listdir(intervention_directory):
#     f = os.path.join(intervention_directory, filename)
#     if os.path.isfile(f):
#         if filename[0:len(output_intervention)] == output_intervention:
#             with open(f, newline='') as csvfile:
#                 strain_counts[intervention][f] = get_strain_counts(list(csv.reader(csvfile)))
# T = max([len(strain_counts[baseline][filename][0]) for filename in strain_counts[baseline]])
# strain_counts_averages = {baseline: np.average([strain_counts[baseline][filename] for filename in strain_counts[baseline]], axis=0),
#                           intervention: np.average([strain_counts[intervention][filename] for filename in strain_counts[intervention]], axis=0)}
# labels = list(range(T))
# strain_1 = 'Recovered'
# strain_2 = 'Recovered'
# plt.figure(figsize=(8, 6))
# font = {'size' : 12}
# plt.rc('font', **font)
# x_lim = 300*144
# plt.xlim(0, x_lim)
# plt.ylim(0, 325000)
# for filename in strain_counts[baseline]:
#     plt.plot(list(range(T)), np.add(strain_counts[baseline][filename][0], strain_counts[baseline][filename][1]), "black", linewidth=1, alpha=0.1)
# plt.plot(list(range(T)), np.add(strain_counts_averages[baseline][0], strain_counts_averages[baseline][1]), "black", linewidth=1, alpha=1.0, label='Recovered: ' + baseline)
# for filename in strain_counts[intervention]:
#     plt.plot(list(range(T)), np.add(strain_counts[intervention][filename][0], strain_counts[intervention][filename][1]), "black", linestyle='dotted', linewidth=1, alpha=0.1)
# plt.plot(list(range(T)), np.add(strain_counts_averages[intervention][0], strain_counts_averages[intervention][1]), "black", linestyle='dotted', linewidth=1, alpha=1.0, label='Recovered: ' + intervention)
# plt.xlabel('Day')
# plt.xticks(ticks=[a*28*144 for a in range((x_lim // (28*144)) + 1)],
#            labels=[a*28 for a in range((x_lim // (28*144)) + 1)])
# plt.legend(loc='upper right')
# plt.grid(False)
# plt.savefig('fig9right.png', bbox_inches='tight')
