from . model import Model

"""
opentiva.alfentanil
===================

This module contains the classes for alfentanil models.

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


class Maitre(Model):
    """Maitre class holds pharmacokinetic parameters for the Maitre alfentanil
    model.

    Reference: PMID: 3099604

    Keo PMID: 1824743 DOI: 10.1097/00000542-199101000-00010
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mcg/ml"
        self.target_unit = "ng/ml"
        self.age_lower = 25
        self.age_upper = 53
        self.weight_lower = 41
        self.weight_upper = 95
        self.pmid = "3099604"
        self.doi = ""
        self.validate_anthropometric_values()

        if sex == 0:
            self.v1 = 0.11 * weight
        elif sex == 1:
            self.v1 = 0.11 * 1.15 * weight

        if age < 40:
            self.k10 = 0.356 / self.v1
        elif age >= 40:
            self.k10 = (0.356 - (0.00269 * (age - 40))) / self.v1

        self.k12 = 0.104
        self.k13 = 0.0170
        self.k21 = 0.0673

        if age < 40:
            self.k31 = 0.0126
        elif age >= 40:
            self.k31 = 0.0126 - (0.000113 * (age - 40))

        self.ke0 = 0.77


class Goresky(Model):
    """Goresky class holds pharmacokinetic parameters for the Goresky
    alfentanil model.

    Reference: PMID: 3118743 DOI: 10.1097/00000542-198711000-00007

    Keo PMID: 1824743 DOI: 10.1097/00000542-199101000-00010
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 2
        self.concentration_unit = "mcg/ml"
        self.target_unit = "ng/ml"
        self.age_lower = 1
        self.age_upper = 14
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "3118743"
        self.doi = "10.1097/00000542-198711000-00007"
        self.validate_anthropometric_values()

        self.v1 = 0.206 * weight
        self.k10 = 0.038
        self.k12 = 0.018
        self.k21 = 0.018

        self.ke0 = 0.77


class Scott(Model):
    """Scott class holds pharmacokinetic parameters for the Scott alfentanil
    model.

    Reference: PMID: 3100765

    Keo PMID: 1824743 DOI: 10.1097/00000542-199101000-00010
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mcg/ml"
        self.target_unit = "ng/ml"
        self.age_lower = -1
        self.age_upper = -1
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "3100765"
        self.doi = ""
        self.validate_anthropometric_values()

        self.v1 = 2.185/70 * weight
        self.cl1 = 0.195

        self.k10 = self.cl1 / self.v1
        self.k12 = 0.656
        self.k13 = 0.113
        self.k21 = 0.214
        self.k31 = 0.017

        self.ke0 = 0.77
