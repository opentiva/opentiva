from . model import Model

"""
opentiva.etomidate
====================

This module contains the classes for the etomidate models.

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


class Kaneda(Model):
    """Kaneda class holds pharmacokinetic parameters for the Kaneda etomidate
    model.

    Reference: PMID: 20498288 DOI: 10.1177/0091270010369242
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 2
        self.concentration_unit = "mg/ml"
        self.age_lower = 18
        self.age_upper = 55
        self.weight_lower = 50
        self.weight_upper = 80
        self.target_unit = "ug/ml"
        self.pmid = "20498288"
        self.doi = "10.1177/0091270010369242"
        self.validate_anthropometric_values()

        self.v1 = 4.45
        self.v2 = 74.9

        self.cl1 = 0.63
        self.cl2 = 3.16

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k21 = self.cl2 / self.v2

        self.ke0 = 0.447


class Lin(Model):
    """Lin class holds pharmacokinetic parameters for the Lin etomidate model.

    Reference: PMID: 21917057 DOI: 10.1111/j.1460-9592.2011.03696.x
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mg/ml"
        self.age_lower = -1
        self.age_upper = 14
        self.weight_lower = -1
        self.weight_upper = -1
        self.target_unit = "ug/ml"
        self.pmid = "21917057"
        self.doi = "10.1111/j.1460-9592.2011.03696.x"
        self.validate_anthropometric_values()

        self.v1 = 9.51 * (age / 4) ** -0.451 * weight / 70
        self.v2 = 11.0 * (weight / 70)
        self.v3 = 79.2 * (age / 4) ** -0.230 * weight / 70

        self.cl1 = 1.50 * (1 - (age - 4) * 0.0288) * (weight / 70) ** 0.75
        self.cl2 = 1.95 * (weight / 70) ** 0.75
        self.cl3 = 1.23 * (weight / 70) ** 0.75

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k13 = self.cl3 / self.v1
        self.k21 = self.cl2 / self.v2
        self.k31 = self.cl3 / self.v3

        self.ke0 = 0.561  # Tpeak 1.5min

