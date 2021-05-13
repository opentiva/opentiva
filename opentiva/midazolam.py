from . model import Model

"""
opentiva.midazolam
===================
This module contains the classes for midazolam models.
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


class Albrecht(Model):
    """Albrecht class holds pharmacokinetic parameters for the Albrecht midazolam
    model.

    Reference: PMID: 10391668 DOI: 10.1016/S0009-9236(99)90084-X

    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = 18
        self.age_upper = -1
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "10391668"
        self.doi = "10.1016/S0009-9236(99)90084-X"
        self.validate_anthropometric_values()

        if age < 65:
            self.v1 = 7.9
            self.cl1 = 399 / 1000

            self.k10 = self.cl1 / self.v1
            self.k12 = 0.19
            self.k13 = 0.065
            self.k21 = 0.06
            self.k31 = 0.0083
            self.ke0 = 0.11
        else:
            self.v1 = 8.5
            self.cl1 = 388 / 1000

            self.k10 = self.cl1 / self.v1
            self.k12 = 0.1
            self.k13 = 0.066
            self.k21 = 0.051
            self.k31 = 0.0069
            self.ke0 = 0.08
