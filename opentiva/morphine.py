from .model import Model

"""
opentiva.morphine
=================

This module contains the classes for morphine models.

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


class Sarton(Model):
    """Sarton class holds pharmacokinetic parameters for the Sarton morphine
    model.

    Reference: PMID: 11046213 DOI: 10.1097/00000542-200011000-00018
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = -1
        self.age_upper = -1
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "11046213"
        self.doi = "10.1097/00000542-200011000-00018"
        self.validate_anthropometric_values()

        self.v1 = 0.075 * weight
        self.k10 = 0.3
        self.k12 = 0.183
        self.k13 = 0.29
        self.k21 = 0.087
        self.k31 = 0.013

        if sex == 0:
            self.ke0 = 0.0073
        else:
            self.ke0 = 0.0024
