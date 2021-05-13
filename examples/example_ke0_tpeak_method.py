#!/usr/bin/env python3
# encoding: utf-8

import opentiva.propofol as propofol


sex = 0  # male
age = 30  # years
weight = 70  # kg
height = 170  # cm

propofol_marsh = propofol.MarshDiprifusor(sex=sex, age=age,
                                          weight=weight,
                                          height=height)

dose = 1 # bolus dose in mg
tpeak = 236 # time to peak effect site concentration in seconds
ce_tpeak = 0.25831 # effect site concentration at tpeak

ke0 = propofol_marsh.ke0_tpeak_method(dose=dose, tpeak=tpeak,
                                      ce_tpeak=ce_tpeak)

print(f"Marsh model Ke0 is {ke0} /second and {ke0/60} /minute")
