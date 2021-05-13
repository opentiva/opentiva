from . model import Model

"""
opentiva.atracurium
===================

This module contains the classes for atracurium models.
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
burns
   bool, true if burn injury present, Marathe only

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


class Fisher(Model):
    """Fisher class holds pharmacokinetic parameters for the Fisher atracurium
    model.

    Reference: PMID: 2360737 DOI: 10.1097/00000542-199007000-00006

    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 2
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = 0
        self.age_upper = 100
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "2360737"
        self.doi = "10.1097/00000542-199007000-00006"
        self.validate_anthropometric_values()

        if age <= 1:
            self.v1 = 100 / 1000 * weight
            self.v2 = 210 / 1000 * weight
            self.cl1 = 3 / 1000 * weight
            self.cl2 = 4.8 / 1000 * weight
            self.k20 = 0.023
            self.ke0 = 0.188
        elif age > 1 and age <= 4:
            self.v1 = 63 / 1000 * weight
            self.v2 = 129 / 1000 * weight
            self.cl1 = 4.2 / 1000 * weight
            self.cl2 = 2.6 / 1000 * weight
            self.k20 = 0.020
            self.ke0 = 0.159
        else:
            self.v1 = 32 / 1000 * weight
            self.v2 = 100 / 1000 * weight
            self.cl1 = 2.8 / 1000 * weight
            self.cl2 = 2.5 / 1000 * weight
            self.k20 = 0.025
            self.ke0 = 0.116

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k21 = self.cl2 / self.v2


class Marathe(Model):
    """Marathe class holds pharmacokinetic parameters for the Marathe atracurium
    model.

    Reference: PMID: 2719307 DOI: 10.1097/00000542-198905000-00007

    """

    def __init__(self, sex: int, age: float, weight: float, height: float,
                 burns: bool = False):
        super().__init__(sex, age, weight, height)

        self.compartments = 1
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = 16
        self.age_upper = 52
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "2719307"
        self.doi = "10.1097/00000542-198905000-00007 "
        self.validate_anthropometric_values()

        if burns:
            self.v1 = 60.9 / 1000 * weight
            self.cl1 = 5.34 / 1000 * weight
            self.ke0 = 0.1
        else:
            self.v1 = 66.3 / 1000 * weight
            self.cl1 = 5.81 / 1000 * weight
            self.ke0 = 0.074

        self.k10 = self.cl1 / self.v1
#       self.k12 = self.cl2 / self.v1
#       self.k21 = self.cl2 / self.v2
