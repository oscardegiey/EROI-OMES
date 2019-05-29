#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    @author: Thibault Martinelle, Adrien Couplet, Oscar de Giey
"""

import Case
import time, os, pstats, cProfile, pickle
import matplotlib.pyplot as plt
from colorama import Style
from linear_optimisation import EROI_optimisation

# === Program Parameters ===
simulation_file = "simulations/RM_FULL.yml" # Path to simulation file
log_level       = 2 # Log level: 0 = very little information, 2 = lots of information
load_results    = 0 # If available, the previously computed results are loaded
save_results    = 1 # Save the computed results
save_json       = 0 # Save case in JSON format

def main():
    # === Program ===
    if log_level >= 1:
        print(Style.BRIGHT + "=== Program Parameters ===" + Style.NORMAL)
        if log_level == 2:
            print("Simulation file: \t{0}\nLog level: \t\t{1}\nLoad results: \t\t{2}\nSave results: \t\t{3}".format(simulation_file, log_level, load_results, save_results))
        print('')

    plt.close('all')

    # Set a profiler to find bottlenecks
    pr = cProfile.Profile()
    pr.enable()
    start = time.time()

    case = Case.Case(simulation_file, log_level)
    if case.optiTD:
        results_filename = os.path.join('results', case.name + '_' + str(case.dt) + '_' + str(case.optiTD) + '_' + str(case.nrTD) + '.pickle')
    else:
        results_filename = os.path.join('results', case.name + '_' + str(case.dt) + '_' + str(case.optiTD) + '.pickle')

    if load_results:
        print('\nLoading data from ' + results_filename)
        with open(results_filename, 'rb') as f:
            [case] = pickle.load(f)
    else:
        EROI_optimisation(case)

    if save_results:
        print('\nSave data into ' + results_filename)
        with open(results_filename, 'wb') as f:
            pickle.dump([case], f)

    if save_json:
        case.save_json()

    total_time = time.time() - start
    # End profiler
    pr.disable()
    st = pstats.Stats(pr)
    st.sort_stats("time")

    # Prints the 1% most time-consuming functions
    if log_level == 2:
        print("\n" + Style.BRIGHT + "=== Execution Profiler ===" + Style.NORMAL)
        st.print_stats(0.01)

    # Write here the plots you want
    timeframe = {}
    timeframe['Year'] = [0,8760]        # Whole year
    timeframe['January'] = [0,743]      # January
    timeframe['Jan_W1'] = [0,167]       # W1 of January
    timeframe['July'] = [4344,5087]     # July
    timeframe['Jul_W1'] = [4344,4511]   # W1 of July
    for d in range(365):
        timeframe['D{:d}'.format(d)] = [d*24,d*24+23]
    timeframe['custom'] = [2900,3100]

    plt.rc('font', family='serif')
    fig, (ax1, ax2) = plt.subplots(1,2)
    case.cells[0].display('complex_balance', timeframe['custom'], resolution="dt", ax=ax1)
    case.cells[0].display('simple_balance', timeframe['custom'], resolution="dt", ax=ax2)
    plt.show()

if __name__ == "__main__":
    main()
