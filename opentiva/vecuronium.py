from . model import Model

"""
opentiva.vecuronium
===================

This module contains the classes for vecuronium models.
All Classes have the same parameters and attributes.

Parameters
----------
sex
    0 for male or 1 for female
age
    in years
weight
    in kg
height
   in cm

Additional
~~~~~~~~~~
temperature
   body temperature in degree Celsius, Caldwell only

Attributes
----------
compartments : int
    number of compartments to model; 1, 2 or 3
concentration_unit : str
    drug concentration unit
target_unit : str
    target concentration unit
age_lower : float
    lower age limit of model; -1 if no limit
age_upper : float
    upper age limit of model; -1 if no limit
weight_lower : float
    lower weight limit of model; -1 if no limit
weight_upper : float
    upper weight limit of model; -1 if no limit
pmid : str
    Pubmed ID of model's reference
doi : str
    Digital Object Identifier (DOI) of model's reference
warning : str
    Warnings relating to non-validated anthropometric values
v1 : float
    volume of central compartment
k10 : float
    equilibrium rate constant from compartment 1 to 0
k12 : float
    equilibrium rate constant from compartment 1 to 2
k13 : float
    equilibrium rate constant from compartment 1 to 3
k21 : float
    equilibrium rate constant from compartment 2 to 1
k31 : float
    equilibrium rate constant from compartment 3 to 1
ke0 : float
    effect compartment equilibrium rate constant
"""


class Caldwell(Model):
    """Caldwell class holds pharmacokinetic parameters for the Caldwell
    vecuronium model.

    Reference: PMID: 10638903 DOI: 10.1097/00000542-200001000-00018

    """

    def __init__(self, sex: int, age: float, weight: float, height: float,
                 temperature: float = 37):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = -1
        self.age_upper = -1
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "10638903"
        self.doi = "10.1097/00000542-200001000-00018"
        self.validate_anthropometric_values()

        self.v1 = 39.8 * weight / 1000
        self.v2 = 67.3 * weight / 1000
        self.v3 = 94.2 * weight / 1000

        self.cl1 = 4.21 * weight / 1000
        self.cl2 = 11.8 * weight / 1000
        self.cl3 = 1.63 * weight / 1000

        if sex == 0:
            self.v1 += self.v1 * 0.225
            self.cl1 -= self.cl1 * 0.22

        if (37 - temperature) > 0:
            self.cl1 -= (37 - temperature) * 0.113

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k13 = self.cl3 / self.v1
        self.k21 = self.cl2 / self.v2
        self.k31 = self.cl3 / self.v3

        self.ke0 = -0.639 + 0.023 * temperature

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError("Expected: float or int")
        if value <= 0:
            raise ValueError("Must be greater than 0")
        self._temperature = value


class Wierda(Model):
    """Wierda class holds pharmacokinetic parameters for the Wierda
    vecuronium model.

    Reference: PMID: 1829656 DOI: 10.1007/BF03007578

    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = 18
        self.age_upper = 60
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "1829656"
        self.doi = "10.1007/BF03007578"
        self.validate_anthropometric_values()

        self.v1 = 0.076 * weight

        self.k10 = 0.09
        self.k12 = 0.22
        self.k13 = 0.023
        self.k21 = 0.23
        self.k31 = 0.008

        self.ke0 = 0.378  # tpeak at 210 seconds
