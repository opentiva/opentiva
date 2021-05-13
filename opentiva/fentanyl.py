from numpy import exp

from . model import Model

"""
opentiva.fentanyl
=================

This module contains the classes for fentanyl models.

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


class Shafer(Model):
    """Shafer class holds pharmacokinetic parameters for the Shafer fentanyl
    model.

    Reference: PMID: 2248388 DOI: 10.1097/00000542-199012000-00005
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mcg/ml"
        self.target_unit = "ng/ml"
        self.age_lower = -1
        self.age_upper = -1
        self.weight_lower = 40
        self.weight_upper = 90
        self.pmid = "2248388"
        self.doi = "10.1097/00000542-199012000-00005"
        self.validate_anthropometric_values()

        self.v1 = 6.09
        self.k10 = 0.0827
        self.k12 = 0.471
        self.k13 = 0.225
        self.k21 = 0.102
        self.k31 = 0.006
        self.ke0 = 0.12


class ShaferW80(Model):
    """ShaferW80 class holds pharmacokinetic parameters for the Shafer fentanyl
    model with weight adjustment if weight > 80.

    Reference: PMID: 2248388 DOI: 10.1097/00000542-199012000-00005

    Adjusted weight if > 80kg PMID: 15329584
    DOI: 10.1097/00000542-200409000-00008
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mcg/ml"
        self.target_unit = "ng/ml"
        self.age_lower = -1
        self.age_upper = -1
        self.weight_lower = 40
        self.weight_upper = -1
        self.pmid = "2248388"
        self.doi = "10.1097/00000542-199012000-00005"
        self.validate_anthropometric_values()

        if age > 80:
            weight_adj = 196.4 * exp(-0.025 * weight) - 53.66
            self.v1 = 2.7 * 0.059 * weight_adj
            self.cl1 = 0.223 + 0.00488 * weight_adj
            self.k10 = self.cl1 / self.v1
        else:
            self.v1 = 6.09
            self.k10 = 0.0827

        self.k12 = 0.471
        self.k13 = 0.225
        self.k21 = 0.102
        self.k31 = 0.006
        self.ke0 = 0.12


class Scott(Model):
    """Scott class holds pharmacokinetic parameters for the Scott fentanyl
    model.

    Reference: PMID: 3100765
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

        self.v1 = 12.7
        self.k10 = 0.0452
        self.k12 = 0.315
        self.k13 = 0.154
        self.k21 = 0.079
        self.k31 = 0.00712
        self.ke0 = 0.12


class Ginsberg(Model):
    """Ginsberg class holds pharmacokinetic parameters for the Ginsberg fentanyl
    model.

    Reference: PMID: 8968173 DOI: 10.1097/00000542-199612000-00007

    Keo PMID: 18270231 DOI: 10.1093/bja/aem408
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 2
        self.concentration_unit = "mcg/ml"
        self.target_unit = "ng/ml"
        self.age_lower = 2
        self.age_upper = 11
        self.weight_lower = 9
        self.weight_upper = 35
        self.pmid = "8968173"
        self.doi = "10.1097/00000542-199612000-00007"
        self.validate_anthropometric_values()

        self.v1 = 0.43 * (weight - 19.8) + 5.8
        self.v2 = 6.2 * (age - 6.4) + 34.4

        self.cl1 = 0.01 * (weight - 19.8) + 0.35
        self.cl2 = 0.82

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k21 = self.cl2 / self.v2

        self.ke0 = 0.28
