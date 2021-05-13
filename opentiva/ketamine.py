from . model import Model

"""
opentiva.ketamine
=================

This module contains the classes for ketamine models.
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


class Clements250(Model):
    """Clements250 class holds pharmacokinetic parameters for the Clements
    ketamine model based on 250 mcg/kg bolus data.

    Reference: PMID: 7459184 DOI: 10.1093/bja/53.1.27
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 2
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = -1
        self.age_upper = -1
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "7459184"
        self.doi = "10.1093/bja/53.1.27"
        self.validate_anthropometric_values()

        self.v1 = 1.7522 * weight

        self.k10 = 0.0109
        self.k12 = 0.0186
        self.k21 = 0.0137

        self.ke0 = 5.2  # from tpeak 1 min


class Domino(Model):
    """Domino class holds pharmacokinetic parameters for the Domino
    ketamine model.

    Reference: PMID: 7198883 DOI:
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
        self.pmid = "7198883"
        self.doi = ""
        self.validate_anthropometric_values()

        self.v1 = 0.063 * weight

        self.k10 = 0.4381
        self.k12 = 0.5921
        self.k13 = 0.5900
        self.k21 = 0.2470
        self.k31 = 0.0146

        self.ke0 = 0.652  # from tpeak 1 min


class Hijazi(Model):
    """Hijazi class holds pharmacokinetic parameters for the Hijazi
    ketamine model.

    Reference: PMID: 12538370 DOI: 10.1093/bja/aeg028
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 2
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = -1
        self.age_upper = 18
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "12538370"
        self.doi = "10.1093/bja/aeg028"
        self.validate_anthropometric_values()

        self.v1 = 1.08 * weight

        self.k10 = 0.0333
        self.k12 = 0.0088
        self.k21 = 0.0030

        self.ke0 = 4.773  # from tpeak 1 min


class Herd(Model):
    """Herd class holds pharmacokinetic parameters for the Herd
    ketamine model.

    Reference: PMID: 17564643 DOI: 10.1111/j.1460-9592.2006.02145.x
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 2
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = 1
        self.age_upper = 14
        self.weight_lower = 10.8
        self.weight_upper = 74.8
        self.pmid = "17564643"
        self.doi = "10.1111/j.1460-9592.2006.02145.x"
        self.validate_anthropometric_values()

        self.v1 = 38.7 * (weight / 70)
        self.v2 = 102 * (weight / 70)

        self.cl1 = 90 / 60 * (weight / 70) ** 0.75
        self.cl2 = 215 / 60 * (weight / 70) ** 0.75

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k21 = self.cl2 / self.v2

        self.ke0 = 2.995  # from tpeak 1 min


class Hornik(Model):
    """Hornik class holds pharmacokinetic parameters for the Hornik
    ketamine model.

    Reference: PMID: 29677389 DOI: 10.1002/jcph.1116
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 2
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = 1
        self.age_upper = 18
        self.weight_lower = 2
        self.weight_upper = 176
        self.pmid = "29677389"
        self.doi = "10.1002/jcph.1116"
        self.validate_anthropometric_values()

        self.v1 = 32.8 * (weight / 70)
        self.v2 = 152 * (weight / 70)

        self.cl1 = 38.9 / 60 * (weight / 70) ** 0.75
        self.cl2 = 54.9 / 60 * (weight / 70) ** 0.75

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k21 = self.cl2 / self.v2

        self.ke0 = 4.212  # from tpeak 1 min
