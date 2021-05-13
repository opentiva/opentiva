from . model import Model

"""
opentiva.cisatracurium
======================

This module contains the classes for cisatracurium models.
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


class Tran(Model):
    """Tran class holds pharmacokinetic parameters for the Tran
    cisatracurium model.

    Reference: PMID: 9806701 DOI: 10.1097/00000539-199811000-00034
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 2
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = 18
        self.age_upper = 65
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "9806701"
        self.doi = "10.1097/00000539-199811000-00034"
        self.validate_anthropometric_values()

        self.v1 = 0.035 * weight

        self.k10 = 0.118
        self.k12 = 0.053
        self.k21 = 0.185
        self.k20 = 0.0237
        self.ke0 = 0.054


class Imbeault(Model):
    """Imbeault class holds pharmacokinetic parameters for the Imbeault
    cisatracurium model.

    Reference: PMID: 16492821 DOI: 10.1213/01.ane.0000195342.29133.ce
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 2
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = 1
        self.age_upper = 6
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "16492821"
        self.doi = "10.1213/01.ane.0000195342.29133.ce"
        self.validate_anthropometric_values()

        self.v1 = 0.087 * weight
        self.k10 = 0.045
        self.k12 = 0.111
        self.k21 = 0.06
        self.k20 = 0.0237
        self.ke0 = 0.115


class Bergeron(Model):
    """Bergeron class holds pharmacokinetic parameters for the Bergeron
    cisatracurium model.

    Reference: PMID: 11506100 DOI: 10.1097/00000542-200108000-00010
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 2
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = 18
        self.age_upper = 65
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "11506100"
        self.doi = "10.1097/00000542-200108000-00010"
        self.validate_anthropometric_values()

        self.v1 = (63.2 + 54.8 + 60.3) / 3 / 1000 * weight
        self.k10 = (0.0448 + 0.0436 + 0.0368) / 3
        self.k12 = (0.1478 + 0.1411 + 0.1387) / 3
        self.k21 = (0.0458 + 0.0417 + 0.0357) / 3
        self.k20 = 0.0237
        self.ke0 = (0.0675 + 0.0568 + 0.0478) / 3
