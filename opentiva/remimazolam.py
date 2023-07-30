from .biometrics import body_mass_index
from .model import Model

"""
opentiva.remimazolam
====================

This module contains the classes for the remimazolam models.

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
asa3
   bool, true if ASA >= 3, Schmith only
asian
   bool, true if asian, Schmith only

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


class Schmith(Model):
    """Schmith class holds pharmacokinetic parameters for the Schmith
    remimazolam model.

    Reference: PMID: 32585566 DOI: 10.1016/j.jclinane.2020.109899
    """

    def __init__(self, sex: int, age: float, weight: float, height: float,
                 asa_3: bool = False, asian: bool = False):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mg/ml"
        self.age_lower = -1
        self.age_upper = -1
        self.weight_lower = -1
        self.weight_upper = -1
        self.target_unit = "ug/ml"
        self.pmid = "32585566"
        self.doi = "10.1016/j.jclinane.2020.109899"
        self.validate_anthropometric_values()

        bmi = body_mass_index(weight, height)

        self.v1 = 2.92 / 70 * weight
        self.v2 = 19.1 / 70 * weight
        if asa_3:
          self.v1 *= 1 - 0.56
          self.v2 *= 1.22

        self.v3 = 9.81 / 70 * weight

        self.cl1 = 61.6 / 70 * weight / 60
        if sex == 1:
          self.cl1 *= 1.11

        self.cl2 = 22.9 / 70 * weight / 60
        self.cl3 = 69.6 / 70 * weight / 60

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k13 = self.cl3 / self.v1
        self.k21 = self.cl2 / self.v2
        self.k31 = self.cl3 / self.v3

        self.ke0 = 8.08 / 60
        if bmi > 25:
          self.ke0 *= 1.17
        if asian:
          self.ke0 *= 1 - 0.48


class Schuttler(Model):
    """Schuttler class holds pharmacokinetic parameters for the Schuttler
    remimazolam model.

    Reference: PMID: 31972655  DOI: 10.1097/ALN.0000000000003103
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mg/ml"
        self.age_lower = -1
        self.age_upper = -1
        self.weight_lower = -1
        self.weight_upper = -1
        self.target_unit = "ug/ml"
        self.pmid = "31972655"
        self.doi = "10.1097/ALN.0000000000003103"
        self.validate_anthropometric_values()

        self.v1 = 4.7 / 75 * weight
        self.v2 = 14.5
        self.v3 = 15.5

        self.cl1 = 1.14
        self.cl2 = 1.04
        self.cl3 = 0.93

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k13 = self.cl3 / self.v1
        self.k21 = self.cl2 / self.v2
        self.k31 = self.cl3 / self.v3

        self.ke0 = 0.27

