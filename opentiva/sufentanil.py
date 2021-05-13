from . model import Model

"""
opentiva.sufentanil
===================

This module contains the classes for the sufentanil models.

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


class Gepts(Model):
    """Gepts class holds pharmacokinetic parameters for the Gepts sufentanil
    model.

    Reference: PMID: 8533912 DOI: 10.1097/00000542-199512000-00010

    Keo: PMID: 1824743 DOI: 10.1097/00000542-199101000-00010
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.drug = "gepts"
        self.concentration_unit = "mcg/ml"
        self.age_lower = 12
        self.age_upper = -1
        self.weight_lower = -1
        self.weight_upper = -1
        self.target_unit = "ng/ml"
        self.pmid = "8533912"
        self.doi = "10.1097/00000542-199512000-00010"
        self.validate_anthropometric_values()

        self.v1 = 14.3
        self.v2 = 63.1
        self.v3 = 261.6

        self.k10 = 0.0645
        self.k12 = 0.1086
        self.k13 = 0.0229
        self.k21 = 0.0245
        self.k31 = 0.0014
        self.ke0 = 0.17559  # calculated with TTPE of 5.6 min


class Greely(Model):
    """Greely class holds pharmacokinetic parameters for the Greely sufentanil
    model.

    Reference: PMID: 2959170

    Keo: PMID: 1824743 DOI: 10.1097/00000542-199101000-00010
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.drug = "gepts"
        self.concentration_unit = "mcg/ml"
        self.age_lower = -1
        self.age_upper = 18
        self.weight_lower = -1
        self.weight_upper = -1
        self.target_unit = "ng/ml"
        self.pmid = "2959170"
        self.doi = ""
        self.validate_anthropometric_values()

        if age < 3:
            self.v1 = 3.09 * weight
            self.k10 = 0.035
            self.k12 = 0.315
            self.k13 = 0.084
            self.k21 = 0.190
            self.k31 = 0.015
        elif age > 13:
            self.v1 = 2.75 * weight
            self.k10 = 0.032
            self.k12 = 0.157
            self.k13 = 0.057
            self.k21 = 0.130
            self.k31 = 0.012
        else:  # age between 3 - 13
            self.v1 = 2.73 * weight
            self.k10 = 0.043
            self.k12 = 0.290
            self.k13 = 0.060
            self.k21 = 0.196
            self.k31 = 0.016

        self.ke0 = 0.227
