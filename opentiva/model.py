#!/usr/bin/env python3
# encoding: utf-8

import warnings

from . helpers import body_mass_index

"""
opentiva.model
==============

This module contains the parent class for the drug models to provide
validation of the instance variables.

Parameters
----------
sex
    0 for male or 1 for female
age
    in year; greater than 0
weight
    in kg; greater than 0
height
    in cm; greater than 0
"""

class Model:

    def __init__(self, sex: int, age: float, weight: float, height: float):
        self.sex = sex
        self.age = age
        self.weight = weight
        self.height = height
        self.bmi = body_mass_index(weight, height)

        self.age_lower = -1
        self.age_upper = -1
        self.weight_lower = -1
        self.weight_upper = -1
        self.bmi_lower = -1
        self.bmi_upper = -1

        self.warning = ""

    @property
    def sex(self):
        return self._sex

    @sex.setter
    def sex(self, value):
        if not isinstance(value, int):
            raise TypeError("Expected: int")
        if value != 0 and value != 1:
            raise ValueError("Must be 0 (male) or 1 (female)")
        self._sex = value

    @property
    def age(self):
        return self._age

    @age.setter
    def age(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError("Expected: float or int")
        if value <= 0:
            raise ValueError("Must be greater than 0")
        self._age = value

    @property
    def weight(self):
        return self._weight

    @weight.setter
    def weight(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError("Expected: float or int")
        if value <= 0:
            raise ValueError("Must be greater than 0")
        self._weight = value

    @property
    def height(self):
        return self._height

    @height.setter
    def height(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError("Expected: float or int")
        if value <= 0:
            raise ValueError("Must be greater than 0")
        self._height = value

    def validate_anthropometric_values(self):
        """ Method validates the anthropometric values are within the range
        specified within the drug model.

        Warnings will be appended to the class's `self.warning` variable and
        to sys.stderr.
        """

        caution = "Proceeding with non-validated anthropometric values " \
                  "may result result in incorrect calculations and " \
                  "malfunction of opentiva."

        if self.age < self.age_lower and self.age_lower != -1:
            description = f"Age {self.age} yrs is below the model's " \
                          f"validated age of {self.age_lower} yrs"
            warnings.warn(f"{description}\n\n{caution}", UserWarning)
            self.warning += (f"{description}\n")

        if self.age > self.age_upper and self.age_upper != -1:
            description = f"Age {self.age} yrs is above the model's " \
                          f"validated age of {self.age_upper} yrs"
            warnings.warn(f"{description}\n\n{caution}", UserWarning)
            self.warning += (f"{description}\n")

        if self.weight < self.weight_lower and self.weight_lower != -1:
            description = f"Weight {self.weight} kg is below the model's " \
                          f"validated weight of {self.weight_lower} kg"
            warnings.warn(f"{description}\n\n{caution}", UserWarning)
            self.warning += (f"{description}\n")

        if self.weight > self.weight_upper and self.weight_upper != -1:
            description = f"Weight {self.weight} kg is above the model's " \
                          f"validated weight of {self.weight_upper} kg"
            warnings.warn(f"{description}\n\n{caution}", UserWarning)
            self.warning += (f"{description}\n")

        if self.bmi < self.bmi_lower and self.bmi_lower != -1:
            description = f"BMI {self.bmi} kg/m^2 is below the model's " \
                          f"validated bmi of {self.bmi_lower} kg/m^2"
            warnings.warn(f"{description}\n\n{caution}", UserWarning)
            self.warning += (f"{description}\n")

        if self.bmi > self.bmi_upper and self.bmi_upper != -1:
            description = f"BMI {self.bmi} kg/m^2 is above the model's " \
                          f"validated bmi of {self.bmi_upper} kg/m^2"
            warnings.warn(f"{description}\n\n{caution}", UserWarning)
            self.warning += (f"{description}\n")
