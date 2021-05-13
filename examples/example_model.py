#!/usr/bin/env python3
# encoding: utf-8

"""
Example target controlled infusion simulation

This script simulates a target controlled infusion model of the Eleveld
propofol model (without opiates)

Documentation of opentiva can be found here:
"""


import matplotlib.pyplot as plt
import numpy as np

import opentiva.propofol as propofol
import opentiva.pump as pump


def plot_concentrations(data):

    fig, ax = plt.subplots()

    ax.plot(data[:,0] / 60, data[:,1], c='b', label="Plasma concentration")
    ax.plot(data[:,0] / 60, data[:,2], c='r', label="Effect site concentration")

    plt.xlabel("Time (minutes)")
    plt.ylabel("Concentration (ug/ml)")
    plt.title("Eleveld without opiates")
    ax.legend()

    plt.show()


if __name__ == "__main__":
    adult35 = propofol.Eleveld(sex=0, age=35, weight=70,
                               height=170, opiates_coadministered=False)

    p1 = pump.Pump(model=adult35,
                   drug_concentration=10,
                   end_time=(1*60*60),
                   maintenance_infusion_duration=10,
                   maintenance_infusion_multiplier=2,
                   cp_limit = 1.2,
                   cp_limit_duration = 10,
                   max_infusion_rate = 1200,
                   bolus_time = 20)

    # Target Ce 3 using revised method at time 0
    p1.add_target(start=0,
                  target=3,
                  duration=10,
                  effect=True,
                  cp_limit=5.0,
                  cp_limit_duration=10,
                  ce_bolus_only=False)

    # Reduce to Ce 2.5 at 5 minutes (300 seconds)
    p1.add_target(start=300,
                  target=2.5,
                  duration=10,
                  effect=True)

    # Target Cp 4 at time 10 minutes (600 seconds) over 1 minute
    p1.add_target(start=1200,
                  target=4,
                  duration=60,
                  effect=False)

    # Reduce to Cp 2 at time 50 minutes (3000 seconds) then stop the
    # maintenance infusions
    p1.add_target(start=3000,
                  target=2,
                  duration=10,
                  effect=False,
                  maintenance_infusions=False)

    concentrations = p1.run()

    plot_concentrations(concentrations)

    # Infusion List
    np.set_printoptions(precision=3)
    np.set_printoptions(suppress=True)

    print("\nInfusion list"
          "\n============="
          "\nColumn 0: Start time (seconds)"
          "\nColumn 1: Dose (mg) per second"
          "\nColumn 2: Duration (seconds)"
          "\nColumn 3: End time (seconds)\n")

    print(p1.infusion_list)

    # Rates list
    np.set_printoptions(precision=1)

    rates = p1.generate_rates_array()

    print("\nRates"
          "\n====="
          "\nColumn 0: Time (seconds)"
          "\nColumn 1: Rate (ml/hr)\n")

    print(rates)
