#!/usr/bin/env python3
# encoding: utf-8

"""
Example user defined infusion

This script simulates a user defined infusion of remifentanil at 10 ml/hr
using the Minto model.

Documentation of opentiva can be found here:
https://opentiva.readthedocs.io
"""


import matplotlib.pyplot as plt

import opentiva.remifentanil as remifentanil
import opentiva.pump as pump


def plot_concentrations(data):

    fig, ax = plt.subplots()

    ax.plot(data[:,0] / 60, data[:,1], c='b', label="Plasma concentration")
    ax.plot(data[:,0] / 60, data[:,2], c='r', label="Effect site concentration")

    plt.xlabel("Time (minutes)")
    plt.ylabel("Concentration (ng/ml)")
    plt.title("Remifentail Minto User Defined Infusion")
    ax.legend()

    plt.show()


if __name__ == "__main__":
    adult75 = remifentanil.Minto(sex=1, age=75, weight=100, height=170)

    p1 = pump.Pump(model=adult75,
                   drug_concentration=50,
                   end_time=(1*60*60),
                   maintenance_infusion_duration=10,
                   maintenance_infusion_multiplier=2,
                   cp_limit = 1.2,
                   cp_limit_duration = 10,
                   max_infusion_rate = 1200,
                   bolus_time = 20)

    total_dose = 50 * 10 * 1 # 50 mcg/ml at 10 ml/hr over 1 hour
    dose_second = total_dose / (1 * 60 * 60)  # Convert to mcg / second

    p1.add_infusion(start=0, dose=dose_second, duration=(1*60*60))

    concentrations = p1.run()

    plot_concentrations(concentrations)

