from collections import defaultdict
import csv
import math
import copy
from numpy.core.fromnumeric import shape
from numpy import genfromtxt
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

def sim(config):
    """Multi-region multi-strain SIR model, solved using the forward Euler method"""

    time_horizon_days = config['time_horizon_days']
    step_size         = config['euler_scheme_step_size_days']
    number_of_regions = config['number_of_regions']
    population_sizes  = config['population_sizes']
    number_of_strains = config['number_of_strains']
    initial_cases     = config['initial_cases']
    beta              = config['transmission_probabilities']
    gamma             = config['removal_rates_per_day']

    assert time_horizon_days     >= 1
    assert number_of_regions     >= 1
    assert len(population_sizes) == number_of_regions
    assert number_of_strains     >= 1
    assert len(initial_cases)    == number_of_regions
    assert initial_cases.shape   == (number_of_regions, number_of_strains)
    assert len(beta)             == number_of_strains
    assert len(gamma)            == number_of_strains

    intra_regional_mixing_per_day     = config['intra_regional_mixing_per_day']
    inter_regional_mixing_per_day     = config['inter_regional_mixing_per_day']
    lockdown                          = config['lockdown']
    lockdown_factor                   = config['lockdown_factor']
    lockdown_start_time_days          = config['lockdown_start_time_days']
    lockdown_end_time_days            = config['lockdown_end_time_days']
    border_closure                    = config['border_closure']
    border_closure_factor             = config['border_closure_factor']
    border_closure_start_time_days    = config['border_closure_start_time_days']
    border_closure_end_time_days      = config['border_closure_end_time_days']
    quarantine                        = config['quarantine']
    quarantine_factor                 = config['quarantine_factor']
    quarantine_factor_start_time_days = config['quarantine_factor_start_time_days']
    quarantine_factor_end_time_days   = config['quarantine_factor_end_time_days']

    time_horizon_ticks                  = int(time_horizon_days * (1/step_size))
    lockdown_start_time_ticks           = int(lockdown_start_time_days * (1/step_size))
    lockdown_end_time_ticks             = int(lockdown_end_time_days * (1/step_size))
    border_closure_start_time_ticks     = int(border_closure_start_time_days * (1/step_size))
    border_closure_end_time_ticks       = int(border_closure_end_time_days * (1/step_size))
    quarantine_factor_start_time_ticks  = int(quarantine_factor_start_time_days * (1/step_size))
    quarantine_factor_end_timee_ticks   = int(quarantine_factor_end_time_days * (1/step_size))

    contact_rate =\
        make_contact_rate_baseline(number_of_regions,
                                   intra_regional_mixing_per_day,
                                   inter_regional_mixing_per_day,
                                   time_horizon_ticks)

    removal_rate =\
        make_removal_rate_baseline(number_of_regions, time_horizon_ticks)

    if lockdown:
        contact_rate =\
            make_contact_rate_lockdown(contact_rate,
                                       lockdown_start_time_ticks,
                                       lockdown_end_time_ticks,
                                       lockdown_factor)
    if border_closure:
        contact_rate =\
            make_contact_rate_border_closure(contact_rate,
                                             border_closure_start_time_ticks,
                                             border_closure_end_time_ticks,
                                             border_closure_factor)

    if quarantine:
        removal_rate =\
            make_removal_rate_quarantine(removal_rate,
                                         quarantine_factor_start_time_ticks,
                                         quarantine_factor_end_timee_ticks,
                                         quarantine_factor)

    # Parameters
    T = time_horizon_ticks
    N = population_sizes
    S = np.zeros((T, number_of_regions), dtype=np.float64)
    I = np.zeros((T, number_of_regions, number_of_strains), dtype=np.float64)
    R = np.zeros((T, number_of_regions), dtype=np.float64)
    g = contact_rate
    h = removal_rate

    # Initial conditions
    for k in range(number_of_regions):
        R[0][k] = 0
        for i in range(number_of_strains):
            I[0][k][i] = initial_cases[k][i]
        S[0][k] = N[k] - sum([I[0][k][i] for i in range(number_of_strains)]) - R[0][k]

    # One-strain dynamics
    # for t in range(T-1):
    #     writer.writerow([I[t][0][0]])
    #     S[t+1][0]    = S[t][0]    - step_size*g[t][0]*beta[0]*I[t][0][0]*S[t][0]/N[0]
    #     I[t+1][0][0] = I[t][0][0] + step_size*g[t][0]*beta[0]*I[t][0][0]*S[t][0]/N[0] \
    #                               - step_size*h[t][0]*gamma[0]*I[t][0][0]
    #     R[t+1][0]    = R[t][0]    + step_size*h[t][0]*gamma[0]*I[t][0][0]

    # Two-strain dynamics (generalizes one-strain dynamics)
    # for t in range(T-1):
    #     S[t+1][0]    = S[t][0]    - step_size*g[t][0][0]*beta[0]*I[t][0][0]*S[t][0]/N[0]\
    #                               - step_size*g[t][0][0]*beta[1]*I[t][0][1]*S[t][0]/N[0]
    #     I[t+1][0][0] = I[t][0][0] + step_size*g[t][0][0]*beta[0]*I[t][0][0]*S[t][0]/N[0]\
    #                               - step_size*h[t][0]*gamma[0]*I[t][0][0]
    #     I[t+1][0][1] = I[t][0][1] + step_size*g[t][0][0]*beta[1]*I[t][0][1]*S[t][0]/N[0]\
    #                               - step_size*h[t][0]*gamma[1]*I[t][0][1]
    #     R[t+1][0]    = R[t][0]    + step_size*h[t][0]*gamma[0]*I[t][0][0]\
    #                               + step_size*h[t][0]*gamma[1]*I[t][0][1]

    # Multi-strain multi-region dynamics (generalizes two-strain dynamics)
    N_bar = np.zeros((number_of_regions, number_of_strains), dtype=float)
    for k in range(number_of_regions):
        for i in range(number_of_strains):
            N_bar[k][i] = 1/N[k]
    for t in range(T-1):
        gtIt = np.matmul(g[t], np.multiply(I[t], N_bar))
        S[t+1] = S[t] - step_size*np.multiply(np.matmul(gtIt, beta), S[t])
        I[t+1] = I[t] + step_size*np.multiply(gtIt, np.tensordot(S[t], beta, axes=0))\
                      - step_size*np.multiply(I[t], np.tensordot(h[t], gamma, axes=0))
        R[t+1] = R[t] + step_size*np.multiply(h[t], np.matmul(I[t], gamma))

    return S, I, R

def make_contact_rate_baseline(number_of_regions, diagonal_mixing, off_diagonal_mixing, T):
    """Makes baseline contact rate matrices with some mixing between regions"""
    basic = np.full((number_of_regions,number_of_regions), off_diagonal_mixing, dtype=np.float64)
    np.fill_diagonal(basic, diagonal_mixing)
    contact_rate_baseline = np.array([basic for t in range(T)], dtype=np.float64)
    return contact_rate_baseline

def make_contact_rate_lockdown(contact_rate_baseline, start_time, end_time, factor):
    """Multiplies contact rates on diagonal"""
    time_horizon = len(contact_rate_baseline)
    modified = np.array([contact_rate_baseline[t] for t in range(time_horizon)], dtype=np.float64)
    shape = contact_rate_baseline[0].shape
    multiplier = np.ones(shape, dtype=np.float64)
    np.fill_diagonal(multiplier, factor)
    for t in range(start_time, end_time):
        modified[t] = np.multiply(contact_rate_baseline[t], multiplier)
    return modified

def make_contact_rate_border_closure(contact_rate_baseline, start_time, end_time, factor):
    """Multiplies contact rates off diagonal"""
    time_horizon = len(contact_rate_baseline)
    modified = np.array([contact_rate_baseline[t] for t in range(time_horizon)], dtype=np.float64)
    shape = contact_rate_baseline[0].shape
    multiplier = np.full(shape, factor, dtype=np.float64)
    np.fill_diagonal(multiplier, 1)
    for t in range(start_time, end_time):
        modified[t] = np.multiply(contact_rate_baseline[t], multiplier)
    return modified

def make_removal_rate_baseline(number_of_regions, T):
    """Makes baseline removal rate arrays for multiple regions"""
    basic = np.ones(number_of_regions)
    removal_rate_baseline = np.array([basic for t in range(T)], dtype=np.float64)
    return removal_rate_baseline

def make_removal_rate_quarantine(removal_rate_baseline, start_time, end_time, factor):
    """Multiplies removal rates"""
    time_horizon = len(removal_rate_baseline)
    modified = np.array([removal_rate_baseline[t] for t in range(time_horizon)], dtype=np.float64)
    shape = removal_rate_baseline[0].shape
    multiplier = np.full(shape, factor, dtype=float)
    for t in range(start_time, end_time):
        modified[t] = np.multiply(removal_rate_baseline[t], multiplier)
    return modified

def calculate_R_differences(alphas, deltas, config_baseline, config_intervention):
    """Simulates baseline and intervention scenarios for variations of beta and gamma and calculates
    final difference in recovered"""

    R_differences = []

    for alpha in tqdm(alphas):
        R_differences_alpha = []
        for delta in deltas:
            # The baseline config
            config_baseline_copy = copy.deepcopy(config_baseline)
            config_baseline_copy['transmission_probabilities'][1] *= alpha
            config_baseline_copy['removal_rates_per_day'][1] *= 1 / delta
            # The baseline scenario
            S_baseline, I_baseline, R_baseline = sim(config_baseline_copy)
            # The intervention config
            config_intervention_copy = copy.deepcopy(config_intervention)
            config_intervention_copy['transmission_probabilities'][1] *= alpha
            config_intervention_copy['removal_rates_per_day'][1] *= 1 / delta
            # The intervention scenario
            S_intervention, I_intervention, R_intervention = sim(config_intervention_copy)
            # Calculate different between the baseline and intervention scenarios
            assert I_baseline.shape == I_intervention.shape
            number_of_regions = I_baseline.shape[1]
            R_differences_alpha.append(sum([R_baseline[-1][k] for k in range(number_of_regions)])\
                                - sum([R_intervention[-1][k] for k in range(number_of_regions)]))
        R_differences.append(R_differences_alpha)

    np.savetxt('output/sir_R_differences.csv', R_differences, delimiter=',')

    return np.array(R_differences, dtype=np.float64)

def plot_R_differences(alphas, deltas, R_differences, show_fig, save_fig):
    """Plots R_differences as a heat map"""

    fig, ax = plt.subplots()
    fig.set_size_inches(16, 16)

    cmap = plt.get_cmap('PuOr') # RdYlGn
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

    R_differences = np.flip(R_differences, axis=0)

    bound = max(abs(np.min(R_differences)), abs(np.max(R_differences)))
    bound = 150000
    bounds = np.arange(-bound, bound, math.ceil((2*bound + 1)/256))
    idx=np.searchsorted(bounds,0)
    bounds=np.insert(bounds,idx,0)
    norm = BoundaryNorm(bounds, cmap.N)

    deltas.reverse()

    ax.set_xticks(np.arange(len(alphas)))
    ax.set_yticks(np.arange(len(deltas)))
    ax.set_xticklabels(alphas, fontsize=18)
    ax.set_yticklabels(deltas, fontsize=18)

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

def sim_and_plot_differences(plot_I, plot_R, config_baseline, config_intervention, show_fig,
                             save_fig, y_lim, intervention_label):
    """Plots output of a simulation"""

    plt.figure(figsize=(8, 6))

    font = {'size' : 12}

    plt.rc('font', **font)

    S_baseline, I_baseline, R_baseline = sim(config_baseline)
    S_intervention, I_intervention, R_intervention = sim(config_intervention)

    assert I_baseline.shape == I_intervention.shape
    assert config_baseline['euler_scheme_step_size_days'] ==\
           config_intervention['euler_scheme_step_size_days']

    T = I_baseline.shape[0]
    step_size = config_baseline['euler_scheme_step_size_days']
    number_of_regions = I_baseline.shape[1]
    number_of_strains = I_baseline.shape[2]

    if number_of_strains == 1:
        strain_colors = ['black']
    else:
        strain_colors = ['red', 'blue', 'green', 'orange']

    if plot_I:
        for k in range(0, number_of_regions):
            for i in range(0, number_of_strains):
                if number_of_strains > 1:
                    plot_label = 'Infected: Baseline: Strain ' + str(i+1)
                else:
                    plot_label = 'Infected: Baseline'
                plt.plot(list(range(T)),
                         [I_baseline[t][k][i] for t in range(T)],
                         strain_colors[i%len(strain_colors)],
                         linewidth=1,
                         alpha=1.0,
                         label=plot_label)
        for k in range(0, number_of_regions):
            for i in range(0, number_of_strains):
                if number_of_strains > 1:
                    plot_label = 'Infected: ' + intervention_label + ': Strain ' + str(i+1)
                else:
                    plot_label = 'Infected: ' + intervention_label
                plt.plot(list(range(T)),
                         [I_intervention[t][k][i] for t in range(T)],
                         strain_colors[i%len(strain_colors)],
                         linewidth=1,
                         linestyle='dotted',
                         alpha=1.0,
                         label=plot_label)

    if plot_R:
        for k in range(0, number_of_regions):
            plt.plot(list(range(T)),
                     [R_baseline[t][k] for t in range(T)],
                     'black',
                     linewidth=1,
                     alpha=1.0,
                     label='Recovered: Baseline')
        for k in range(0, number_of_regions):
            plt.plot(list(range(T)),
                     [R_intervention[t][k] for t in range(T)],
                     'black',
                     linewidth=1,
                     linestyle='dotted',
                     alpha=1.0,
                     label='Recovered: ' + intervention_label)

    plt.xlabel('Day')
    plt.xticks(ticks=[a*int(28/step_size) for a in range((T // int(28/step_size)) + 1)],
               labels=[a*int(28) for a in range((T // int(28/step_size)) + 1)])
    plt.grid(False)
    plt.xlim([0, T])
    if y_lim is not None:
        plt.ylim([0, y_lim])

    plt.legend(loc='upper right')

    if show_fig:
        plt.show()
    if save_fig:
        plt.savefig('sim_output.png', bbox_inches='tight')

def sim_and_plot_lockdown_start_times(alpha, delta, config_baseline_2_strain,
                                                    config_lockdown_2_strain, final_start_day):
    """Varies lockdown start time and plots differences"""
    R_differences = []
    config_baseline_2_strain['transmission_probabilities'][1] *= alpha
    config_baseline_2_strain['removal_rates_per_day'][1] *= 1 / delta
    config_lockdown_2_strain['transmission_probabilities'][1] *= alpha
    config_lockdown_2_strain['removal_rates_per_day'][1] *= 1 / delta
    for i in tqdm(range(final_start_day)):
        config_baseline_2_strain['lockdown_start_time_days'] = i
        config_baseline_2_strain['lockdown_end_time_days'] = i + 28
        config_lockdown_2_strain['lockdown_start_time_days'] = i
        config_lockdown_2_strain['lockdown_end_time_days'] = i + 28
        S_baseline, I_baseline, R_baseline = sim(config_baseline_2_strain)
        S_intervention, I_intervention, R_intervention = sim(config_lockdown_2_strain)
        assert I_baseline.shape == I_intervention.shape
        number_of_regions = I_baseline.shape[1]
        R_differences.append(sum([R_baseline[-1][k] for k in range(number_of_regions)])\
                           - sum([R_intervention[-1][k] for k in range(number_of_regions)]))

    plt.figure(figsize=(8, 6))

    font = {'size' : 12}

    plt.rc('font', **font)

    plot_label = 'Reduction in Total Cumulative Cases'
    plt.plot(list(range(final_start_day)),
             [R_differences[t] for t in range(final_start_day)],'black',
             linewidth=1, alpha=1.0, label=plot_label)
    plt.xlabel('Start day of 28 day lockdown')
    plt.xticks(ticks=[a*int(7) for a in range((final_start_day // int(7)) + 1)],
               labels=[a*int(7) for a in range((final_start_day // int(7)) + 1)])
    plt.grid(False)
    plt.xlim([0, final_start_day])
    # plt.ylim([-140000, 140000])
    plt.legend(loc='upper right')
    plt.savefig('start_time_differences.png', bbox_inches='tight')

config_baseline_1_strain = \
       {'time_horizon_days': 400,
        'euler_scheme_step_size_days': 1/144,
        'number_of_regions': 1,
        'population_sizes': np.array([625960]),
        'number_of_strains': 1,
        'initial_cases': np.array([[320]]),
        'transmission_probabilities': np.array([0.00035]),
        'removal_rates_per_day': np.array([1 / 9]),
        'intra_regional_mixing_per_day': 778,
        'inter_regional_mixing_per_day': 0.0,
        'lockdown': False,
        'lockdown_factor': 0.1,
        'lockdown_start_time_days': 21,
        'lockdown_end_time_days': 49,
        'border_closure': False,
        'border_closure_factor': 0.1,
        'border_closure_start_time_days': 21,
        'border_closure_end_time_days': 49,
        'quarantine': False,
        'quarantine_factor': 10,
        'quarantine_factor_start_time_days': 21,
        'quarantine_factor_end_time_days': 49}

config_lockdown_1_strain =\
       {'time_horizon_days': 400,
        'euler_scheme_step_size_days': 1/144,
        'number_of_regions': 1,
        'population_sizes': np.array([625960]),
        'number_of_strains': 1,
        'initial_cases': np.array([[320]]),
        'transmission_probabilities': np.array([0.00035]),
        'removal_rates_per_day': np.array([1 / 9]),
        'intra_regional_mixing_per_day': 778,
        'inter_regional_mixing_per_day': 0.0,
        'lockdown': True,
        'lockdown_factor': 0.1,
        'lockdown_start_time_days': 21,
        'lockdown_end_time_days': 49,
        'border_closure': False,
        'border_closure_factor': 0.1,
        'border_closure_start_time_days': 21,
        'border_closure_end_time_days': 49,
        'quarantine': False,
        'quarantine_factor': 10,
        'quarantine_factor_start_time_days': 21,
        'quarantine_factor_end_time_days': 49}

config_baseline_2_strain =\
       {'time_horizon_days': 400,
        'euler_scheme_step_size_days': 1/144,
        'number_of_regions': 1,
        'population_sizes': np.array([625960]),
        'number_of_strains': 2,
        'initial_cases': np.array([[320,320]]),
        'transmission_probabilities': np.array([0.00035, 0.00035]),
        'removal_rates_per_day': np.array([1 / 9, 1 / 9]),
        'intra_regional_mixing_per_day': 778,
        'inter_regional_mixing_per_day': 0.0,
        'lockdown': False,
        'lockdown_factor': 0.1,
        'lockdown_start_time_days': 21,
        'lockdown_end_time_days': 49,
        'border_closure': False,
        'border_closure_factor': 0.1,
        'border_closure_start_time_days': 21,
        'border_closure_end_time_days': 49,
        'quarantine': False,
        'quarantine_factor': 10,
        'quarantine_factor_start_time_days': 21,
        'quarantine_factor_end_time_days': 49}

config_lockdown_2_strain =\
       {'time_horizon_days': 400,
        'euler_scheme_step_size_days': 1/144,
        'number_of_regions': 1,
        'population_sizes': np.array([625960]),
        'number_of_strains': 2,
        'initial_cases': np.array([[320,320]]),
        'transmission_probabilities': np.array([0.00035, 0.00035]),
        'removal_rates_per_day': np.array([1 / 9, 1 / 9]),
        'intra_regional_mixing_per_day': 778,
        'inter_regional_mixing_per_day': 0.0,
        'lockdown': True,
        'lockdown_factor': 0.1,
        'lockdown_start_time_days': 21,
        'lockdown_end_time_days': 49,
        'border_closure': False,
        'border_closure_factor': 0.1,
        'border_closure_start_time_days': 21,
        'border_closure_end_time_days': 49,
        'quarantine': False,
        'quarantine_factor': 10,
        'quarantine_factor_start_time_days': 21,
        'quarantine_factor_end_time_days': 49}

config_baseline_2_strain_2_region =\
       {'time_horizon_days': 400,
        'euler_scheme_step_size_days': 1/144,
        'number_of_regions': 2,
        'population_sizes': np.array([625960, 625960]),
        'number_of_strains': 2,
        'initial_cases': np.array([[320,320],[0,0]]),
        'transmission_probabilities': np.array([0.00035, 0.00035]),
        'removal_rates_per_day': np.array([1 / 9, 1 / 9]),
        'intra_regional_mixing_per_day': 778,
        'inter_regional_mixing_per_day': 0.00778,
        'lockdown': False,
        'lockdown_factor': 0.1,
        'lockdown_start_time_days': 21,
        'lockdown_end_time_days': 49,
        'border_closure': False,
        'border_closure_factor': 0.1,
        'border_closure_start_time_days': 21,
        'border_closure_end_time_days': 49,
        'quarantine': False,
        'quarantine_factor': 10,
        'quarantine_factor_start_time_days': 21,
        'quarantine_factor_end_time_days': 49}

config_lockdown_2_strain_2_region =\
       {'time_horizon_days': 400,
        'euler_scheme_step_size_days': 1/144,
        'number_of_regions': 2,
        'population_sizes': np.array([625960, 625960]),
        'number_of_strains': 2,
        'initial_cases': np.array([[320,320],[0,0]]),
        'transmission_probabilities': np.array([0.00035, 0.00035]),
        'removal_rates_per_day': np.array([1 / 9, 1 / 9]),
        'intra_regional_mixing_per_day': 778,
        'inter_regional_mixing_per_day': 0.00778,
        'lockdown': True,
        'lockdown_factor': 0.1,
        'lockdown_start_time_days': 21,
        'lockdown_end_time_days': 49,
        'border_closure': False,
        'border_closure_factor': 0.1,
        'border_closure_start_time_days': 21,
        'border_closure_end_time_days': 49,
        'quarantine': False,
        'quarantine_factor': 10,
        'quarantine_factor_start_time_days': 21,
        'quarantine_factor_end_time_days': 49}

# Figure 1 (left)
# config_baseline_1_strain['time_horizon_days'] = 200
# config_lockdown_1_strain['time_horizon_days'] = 200
# sim_and_plot_differences(True, False, config_baseline_1_strain, config_lockdown_1_strain,
#                          False, True, 150000, 'Lockdown')

# Figure 1 (right)
# config_baseline_1_strain['time_horizon_days'] = 200
# config_lockdown_1_strain['time_horizon_days'] = 200
# sim_and_plot_differences(False, True, config_baseline_1_strain, config_lockdown_1_strain,
#                          False, True, 690000, 'Lockdown')

# Figure 2
# alphas = [round(0.1 + (0.1 * r), 6) for r in range(49)]
# deltas = [round(0.1 + (0.1 * r), 6) for r in range(49)]
# R_differences = calculate_R_differences(alphas, deltas, config_baseline_2_strain,
#                                                         config_lockdown_2_strain)
# R_differences = np.genfromtxt("output/sir_R_differences.csv", delimiter=",", dtype=float)
# plot_R_differences(alphas, deltas, R_differences, False, True)

# Figure 3 (left)
# alpha = 3.0
# delta = 0.2
# config_baseline_2_strain['time_horizon_days'] = 200
# config_lockdown_2_strain['time_horizon_days'] = 200
# config_baseline_2_strain['transmission_probabilities'][1] *= alpha
# config_baseline_2_strain['removal_rates_per_day'][1] *= 1 / delta
# config_lockdown_2_strain['transmission_probabilities'][1] *= alpha
# config_lockdown_2_strain['removal_rates_per_day'][1] *= 1 / delta
# sim_and_plot_differences(True, False, config_baseline_2_strain, config_lockdown_2_strain,
#                          False, True, 100000, 'Lockdown')

# Figure 3 (right)
# alpha = 3.0
# delta = 0.2
# config_baseline_2_strain['time_horizon_days'] = 200
# config_lockdown_2_strain['time_horizon_days'] = 200
# config_baseline_2_strain['transmission_probabilities'][1] *= alpha
# config_baseline_2_strain['removal_rates_per_day'][1] *= 1 / delta
# config_lockdown_2_strain['transmission_probabilities'][1] *= alpha
# config_lockdown_2_strain['removal_rates_per_day'][1] *= 1 / delta
# sim_and_plot_differences(False, True, config_baseline_2_strain, config_lockdown_2_strain,
#                          False, True, 690000, 'Lockdown')

# Figure 4
# alpha = 3.6
# delta = 0.2
# config_baseline_2_strain['initial_cases'] = np.array([[1,1]])
# config_lockdown_2_strain['initial_cases'] = np.array([[1,1]])
# config_baseline_2_strain['removal_rates_per_day'] = np.array([1 / 6, 1 / 6])
# config_lockdown_2_strain['removal_rates_per_day'] = np.array([1 / 6, 1 / 6])
# config_baseline_2_strain['intra_regional_mixing_per_day'] = 1000
# config_lockdown_2_strain['intra_regional_mixing_per_day'] = 1000
# sim_and_plot_lockdown_start_times(alpha, delta, config_baseline_2_strain,
#                                                 config_lockdown_2_strain, 84)

# Figure 5 (left)
# alpha = 3.6
# delta = 0.2
# config_baseline_2_strain['initial_cases'] = np.array([[1,1]])
# config_lockdown_2_strain['initial_cases'] = np.array([[1,1]])
# config_baseline_2_strain['removal_rates_per_day'] = np.array([1 / 6, 1 / 6])
# config_lockdown_2_strain['removal_rates_per_day'] = np.array([1 / 6, 1 / 6])
# config_baseline_2_strain['transmission_probabilities'][1] *= alpha
# config_baseline_2_strain['removal_rates_per_day'][1] *= 1 / delta
# config_lockdown_2_strain['transmission_probabilities'][1] *= alpha
# config_lockdown_2_strain['removal_rates_per_day'][1] *= 1 / delta
# config_baseline_2_strain['intra_regional_mixing_per_day'] = 1000
# config_lockdown_2_strain['intra_regional_mixing_per_day'] = 1000
# config_baseline_2_strain['time_horizon_days'] = 240
# config_lockdown_2_strain['time_horizon_days'] = 240
# config_lockdown_2_strain['lockdown_start_time_days'] = 25
# config_lockdown_2_strain['lockdown_end_time_days'] = 25 + 28
# sim_and_plot_differences(True, False, config_baseline_2_strain, config_lockdown_2_strain,
#                          False, True, 60000, 'Lockdown')

# Figure 5 (right)
# alpha = 3.6
# delta = 0.2
# config_baseline_2_strain['initial_cases'] = np.array([[1,1]])
# config_lockdown_2_strain['initial_cases'] = np.array([[1,1]])
# config_baseline_2_strain['removal_rates_per_day'] = np.array([1 / 6, 1 / 6])
# config_lockdown_2_strain['removal_rates_per_day'] = np.array([1 / 6, 1 / 6])
# config_baseline_2_strain['transmission_probabilities'][1] *= alpha
# config_baseline_2_strain['removal_rates_per_day'][1] *= 1 / delta
# config_lockdown_2_strain['transmission_probabilities'][1] *= alpha
# config_lockdown_2_strain['removal_rates_per_day'][1] *= 1 / delta
# config_baseline_2_strain['intra_regional_mixing_per_day'] = 1000
# config_lockdown_2_strain['intra_regional_mixing_per_day'] = 1000
# config_baseline_2_strain['time_horizon_days'] = 240
# config_lockdown_2_strain['time_horizon_days'] = 240
# config_lockdown_2_strain['lockdown_start_time_days'] = 33
# config_lockdown_2_strain['lockdown_end_time_days'] = 33 + 28
# sim_and_plot_differences(True, False, config_baseline_2_strain, config_lockdown_2_strain,
#                          False, True, 60000, 'Lockdown')

# Figure 6
# alpha = 3.6
# delta = 0.2
# config_baseline_2_strain_2_region['initial_cases'] = np.array([[1,1],[0,0]])
# config_lockdown_2_strain_2_region['initial_cases'] = np.array([[1,1],[0,0]])
# config_baseline_2_strain_2_region['removal_rates_per_day'] = np.array([1 / 6, 1 / 6])
# config_lockdown_2_strain_2_region['removal_rates_per_day'] = np.array([1 / 6, 1 / 6])
# config_baseline_2_strain_2_region['intra_regional_mixing_per_day'] = 1000
# config_lockdown_2_strain_2_region['intra_regional_mixing_per_day'] = 1000
# config_baseline_2_strain_2_region['inter_regional_mixing_per_day'] = 0.001
# config_lockdown_2_strain_2_region['inter_regional_mixing_per_day'] = 0.001
# sim_and_plot_lockdown_start_times(alpha, delta, config_baseline_2_strain_2_region,
#                                                 config_lockdown_2_strain_2_region, 112)
