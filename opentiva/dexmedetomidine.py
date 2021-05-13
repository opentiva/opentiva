from . model import Model

"""
opentiva.dexmedetomidine
========================
This module contains the classes for dexmedetomidine models.
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


class Hannivoort(Model):
    """Hannivoort class holds pharmacokinetic parameters for the Hannivoort
    dexmedetomidine model.

    Reference: PMID: 26068206 DOI: 10.1097/ALN.0000000000000740

    Keo for sedation: PMID: 28854538 DOI: 10.1093/bja/aex085
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mcg/ml"
        self.target_unit = "ng/ml"
        self.age_lower = 20
        self.age_upper = 70
        self.weight_lower = 51
        self.weight_upper = 110
        self.pmid = "26068206"
        self.doi = "10.1097/ALN.0000000000000740"
        self.validate_anthropometric_values()

        self.v1 = 1.78 * (weight / 70)
        self.v2 = 30.3 * (weight / 70)
        self.v3 = 52.0 * (weight / 70)

        self.cl1 = 0.686 * (weight / 70) ** 0.75
        self.cl2 = 2.98 * (self.v2 / 30.3) ** 0.75
        self.cl3 = 0.602 * (self.v3 / 52.0) ** 0.75

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k13 = self.cl3 / self.v1
        self.k21 = self.cl2 / self.v2
        self.k31 = self.cl3 / self.v3

        self.ke0 = 0.0428


class PerezGuille(Model):
    """PerezGuille class holds pharmacokinetic parameters for the Pérez-Guillé
    dexmedetomidine model.

    PMID: 29782406 DOI: 10.1213/ANE.0000000000003413

    Keo for sedation: PMID: 28854538 DOI: 10.1093/bja/aex085
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 2
        self.concentration_unit = "mcg/ml"
        self.target_unit = "ng/ml"
        self.age_lower = 2
        self.age_upper = 18
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "29782406"
        self.doi = "10.1213/ANE.0000000000003413"
        self.validate_anthropometric_values()

        theta_cl1 = 20.8
        theta_cl2 = 75.8
        theta_v1 = 21.9
        theta_v2 = 81.2

        self.v1 = theta_v1 * (weight / 70) ** 0.75
        self.v2 = theta_v2 * (weight / 70) ** 0.75

        self.cl1 = theta_cl1 * (weight / 70) ** 0.75
        self.cl1 /= 60  # convert to L/min
        self.cl2 = theta_cl2 * (weight / 70) ** 0.75
        self.cl2 /= 60  # convert to L/min

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k21 = self.cl2 / self.v2

        self.ke0 = 0.0428


class Rolle(Model):
    """Rolle class holds pharmacokinetic parameters for the Rolle
    dexmedetomidine model.

    Reference: PMID: 29661414 DOI: 10.1016/j.bja.2018.01.040

    Keo for sedation: PMID: 28854538 DOI: 10.1093/bja/aex085
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 2
        self.concentration_unit = "mcg/ml"
        self.target_unit = "ng/ml"
        self.age_lower = 23
        self.age_upper = 59
        self.weight_lower = 47
        self.weight_upper = 126
        self.pmid = "29661414"
        self.doi = "10.1016/j.bja.2018.01.040"
        self.validate_anthropometric_values()

        if sex == 0:
            whs_max = 42.92
            whs_50 = 30.93
        else:
            whs_max = 37.99
            whs_50 = 35.98

        ffm = whs_max * (height / 100) ** 2 * \
              (weight / (whs_50 * (height / 100) ** 2 + weight))

        theta_1 = 30.3
        theta_2 = 71.3
        theta_4 = 0.585
        theta_5 = 1.96

        self.v1 = theta_1 * ffm / 45
        self.v2 = theta_2 * ffm / 45

        self.cl1 = theta_4 * ffm / 45
        self.cl2 = theta_5 * ffm / 45

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k21 = self.cl2 / self.v2

        self.ke0 = 0.0428


class Dyck(Model):
    """Dyck class holds pharmacokinetic parameters for the Dyck
    dexmedetomidine model.

    Reference: PMID: 8098191 DOI: 10.1097/00000542-199305000-00003

    Keo for sedation: PMID: 28854538 DOI: 10.1093/bja/aex085
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mcg/ml"
        self.target_unit = "ng/ml"
        self.age_lower = -1
        self.age_upper = -1
        self.weight_lower = 60
        self.weight_upper = 98
        self.pmid = "8098191"
        self.doi = "10.1093/bja/aex085"
        self.validate_anthropometric_values()

        self.v1 = 7.99
        self.v2 = 13.8
        self.v3 = 187

        self.cl1 = 0.00791 * height - 0.927
        self.cl2 = 2.26
        self.cl3 = 1.99

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k13 = self.cl3 / self.v1
        self.k21 = self.cl2 / self.v2
        self.k31 = self.cl3 / self.v3

        self.ke0 = 0.0428
