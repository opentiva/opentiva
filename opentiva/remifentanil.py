from numpy import exp

from . helpers import lbm_dubois, ffm_alsallami, ffm_janmahasation, \
                      body_mass_index
from . model import Model

"""
opentiva.remifentanil
=====================

This module contains the classes for the remifentanil models.

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


class Minto(Model):
    """Minto class holds the pharmacokinetic parameters for the Minto
    remifentanil model.

    Reference: PMID: 9009935 DOI: 10.1097/00000542-199701000-00004
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mcg/ml"
        self.target_unit = "ng/ml"
        self.age_lower = 12
        self.age_upper = -1
        self.weight_lower = 30
        self.weight_upper = -1
        self.pmid = "9009935"
        self.doi = "10.1097/00000542-199701000-00004"
        self.validate_anthropometric_values()

        lean_body_mass = lbm_dubois(sex, weight, height)

        self.v1 = 5.1 - 0.0201 * (age - 40) + 0.072 * (lean_body_mass - 55)
        self.v2 = 9.82 - 0.0811 * (age - 40) + 0.108 * (lean_body_mass - 55)
        self.v3 = 5.42

        self.k10 = (2.6 - 0.0162 * (age - 40) + 0.0191 *
                    (lean_body_mass - 55)) / self.v1
        self.k12 = (2.05 - 0.0301 * (age - 40)) / self.v1
        self.k13 = (0.076 - 0.00113 * (age - 40)) / self.v1
        self.k21 = self.k12 * (self.v1 / self.v2)
        self.k31 = self.k13 * (self.v1 / self.v3)
        self.ke0 = 0.595 - 0.007 * (age - 40)


class Eleveld(Model):
    """Eleveld class holds the pharmacokinetic parameters for the Eleveld
    remifentanil model.

    Reference: PMID: 28509794 DOI: 10.1097/ALN.0000000000001634
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
        self.pmid = "28509794"
        self.doi = " 10.1097/ALN.0000000000001634"
        self.validate_anthropometric_values()

        # Reference person
        # Age 35, Weight 70kg, Height 170

        v1_ref = 5.81
        v2_ref = 8.82
        v3_ref = 5.03
        cl1_ref = 2.58
        cl2_ref = 1.72
        cl3_ref = 0.124
        theta_1 = 2.88
        theta_2 = -0.00554
        theta_3 = -0.00327
        theta_4 = -0.0315
        theta_5 = 0.470
        theta_6 = -0.0260

        def ageing(x, age):
            return exp(x * (age - 35))

        def sigmoid(x, e50, y):
            return (x ** y) / ((x ** y) + (e50 ** y))

        kmat = sigmoid(weight, theta_1, 2)
        kmat_ref = sigmoid(70, theta_1, 2)

        size = ffm_alsallami(sex, age, weight, height) / \
               ffm_alsallami(0, 35, 70, 170)

        if sex:
            ksex = 1
        else:
            ksex = 1 + theta_5 * sigmoid(age, 12, 6) * \
                   (1 - sigmoid(age, 45, 6))

        self.v1 = v1_ref * size * ageing(theta_2, age)
        self.v2 = v2_ref * size * ageing(theta_3, age) * ksex
        self.v3 = v3_ref * size * ageing(theta_4, age) * \
                  exp(theta_6 * (weight - 70))

        self.cl1 = cl1_ref * size ** 0.75 * (kmat / kmat_ref) * \
                   ksex * ageing(theta_3, age)
        self.cl2 = cl2_ref * (self.v2 / v2_ref) ** 0.75 * \
                   ageing(theta_2, age) * ksex
        self.cl3 = cl3_ref * (self.v3 / v3_ref) ** 0.75 * \
                   ageing(theta_2, age)

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k13 = self.cl3 / self.v1
        self.k21 = self.cl2 / self.v2
        self.k31 = self.cl3 / self.v3
        self.ke0 = 1.09 * ageing(-0.0289, age)


class RigbyJones(Model):
    """RigbyJones class holds the pharmacokinetic parameters for the
    Rigby-Jones remifentanil model.

    Reference: PMID: 17578905 DOI: 10.1093/bja/aem135

    Keo PMID: 18270231 DOI: 10.1093/bja/aem408
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 2
        self.concentration_unit = "mcg/ml"
        self.target_unit = "ng/ml"
        self.age_lower = 1
        self.age_upper = 9
        self.weight_lower = 3
        self.weight_upper = 40
        self.pmid = "17578905"
        self.doi = "10.1093/bja/aem135"
        self.validate_anthropometric_values()

        self.v1 = 716 * (weight / 10.5) ** 0.75 / 1000
        self.v2 = 840 * (weight / 10.5) ** 0.75 / 1000

        self.cl1 = 963 * (weight / 10.5) / 1000
        self.cl2 = 1480 * (weight / 10.5) / 1000

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k21 = self.cl2 / self.v2
        self.ke0 = 0.71


class Kim(Model):
    """Kim class holds the pharmacokinetic parameters for the
    Kim remifentanil model.

    Reference: PMID: 28509796 DOI: 10.1097/ALN.0000000000001635

    Keo PMID: 18270231 DOI: 10.1093/bja/aem408
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mcg/ml"
        self.target_unit = "ng/ml"
        self.age_lower = 20
        self.age_upper = 85
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "28509796"
        self.doi = "10.1097/ALN.0000000000001635"
        self.validate_anthropometric_values()

        def ageing(x, age):
            return exp(x * (age - 35))

        bmi = body_mass_index(weight, height)
        ffm = ffm_janmahasation(sex, weight, height)

        theta_1 = 4.76
        theta_2 = 8.4
        theta_3 = 4
        theta_4 = 2.77
        theta_5 = 1.94
        theta_6 = 0.197
        theta_9 = 0.658
        theta_10 = 0.573
        theta_11 = 0.0936
        theta_12 = 0.0477
        theta_13 = 0.336
        theta_14 = 0.0149
        theta_15 = 0.0280

        self.v1 = theta_1 * (weight / 74.5) ** theta_9
        self.v2 = theta_2 * (ffm / 52.3) ** theta_10 - theta_11
        self.v3 = theta_3 - theta_12 * (age - 37)

        self.cl1 = theta_4 * (weight / 74.5) ** theta_13 - theta_14 * \
                   (age - 37)
        self.cl2 = theta_5 - theta_15 * (age - 37)
        self.cl3 = theta_6

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k13 = self.cl3 / self.v1
        self.k21 = self.cl2 / self.v2
        self.k31 = self.cl3 / self.v3
        self.ke0 = 1.09 * ageing(-0.0289, age)  # Keo from Eleveld model
