from numpy import exp

from . helpers import lbm_dubois, ffm_alsallami
from . model import Model

"""
opentiva.propofol
=================

This module contains the classes for propofol models.

All Classes have the same following parameters and attributes.

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
opiates_coadministered
   bool; true if opiates co-administered, Eleveld model only
bolus_data
   bool; true for bolus data, Schuttler model only
venous_data
   bool; true for venous data, Schuttler model only

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


class MarshDiprifusor(Model):
    """MarshDiprifusor class holds pharmacokinetic parameters for the
       Diprifusor Marsh propofol model with Keo 0.26.

       Reference: PMID: 1859758 DOI: 10.1093/bja/67.1.41
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = 16
        self.age_upper = -1
        self.weight_lower = -1
        self.weight_upper = 150
        self.pmid = "1859758"
        self.doi = "10.1093/bja/67.1.41"
        self.validate_anthropometric_values()

        self.v1 = 0.228 * weight
        self.v2 = 0.463 * weight
        self.v3 = 2.893 * weight

        self.k10 = 0.119
        self.k12 = 0.112
        self.k13 = 0.0419
        self.k21 = 0.055
        self.k31 = 0.0033
        self.ke0 = 0.26


class MarshModified(Model):
    """MarshModified class holds pharmacokinetic parameters for the Modified
       Marsh propofol model with  Keo 1.2.

       Reference: PMID: 1859758 DOI: 10.1093/bja/67.1.41
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = 16
        self.age_upper = -1
        self.weight_lower = -1
        self.weight_upper = 150
        self.pmid = "1859758"
        self.doi = "10.1093/bja/67.1.41"
        self.validate_anthropometric_values()

        self.v1 = 0.228 * weight
        self.v2 = 0.463 * weight
        self.v3 = 2.893 * weight

        self.k10 = 0.119
        self.k12 = 0.112
        self.k13 = 0.0419
        self.k21 = 0.055
        self.k31 = 0.0033
        self.ke0 = 1.2


class Schnider(Model):
    """Schnider class holds pharmacokinetic parameters for the Schnider propofol
    model.

    Reference: PMID: 9605675 DOI: 10.1097/00000542-199805000-00006
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
        if sex == 0:
            self.bmi_upper = 42
        elif sex == 1:
            self.bmi_upper = 35
        self.pmid = "9605675"
        self.doi = "10.1097/00000542-199805000-00006"
        self.validate_anthropometric_values()

        lbm = lbm_dubois(sex, weight, height)

        self.v1 = 4.27
        self.v2 = 18.9 - 0.391 * (age - 52)
        self.v3 = 238

        cl1 = 1.89 + 0.0456 * (weight - 77) - 0.0681 * (lbm - 59) \
              + 0.0264 * (height - 177)
        cl2 = 1.29 - 0.024 * (age - 53)
        cl3 = 0.836

        self.k10 = cl1 / self.v1
        self.k12 = cl2 / self.v1
        self.k13 = cl3 / self.v1
        self.k21 = cl2 / self.v2
        self.k31 = cl3 / self.v3
        self.ke0 = 0.456  # TTPE 1.6 minutes is used in original model


class Paedfusor(Model):
    """Paedfusor class holds pharmacokinetic parameters for the Paedfusor propofol
    model.

    Reference: PMID: 15941735 DOI: 10.1093/bja/aei567
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = 1
        self.age_upper = 16
        self.weight_lower = 5
        self.weight_upper = 61
        self.pmid = "15941735"
        self.doi = "10.1093/bja/aei567"
        self.validate_anthropometric_values()

        self.k12 = 0.114
        self.k13 = 0.0419
        self.k21 = 0.055
        self.k31 = 0.0033
        self.ke0 = 0.26

        if age <= 12:
            self.v1 = 458.4 * weight / 1000
            self.k10 = 0.1527 * (weight ** -0.3)
        elif age >= 13:
            self.v1 = 400 * weight / 1000
            self.k10 = 0.0678
        elif age >= 14:
            self.v1 = 342 * weight / 1000
            self.k10 = 0.0792
        elif age >= 15:
            self.v1 = 284 * weight / 1000
            self.k10 = 0.0954
        elif age >= 16:
            self.v1 = 228.57 * weight / 1000
            self.k10 = 0.119

        self.v2 = self.v1 * self.k12 / self.k21 / 1000
        self.v3 = self.v1 * self.k13 / self.k31 / 1000


class Kataria(Model):
    """Kataria class holds pharmacokinetic parameters for the Kataria propofol
    model.

    Reference: PMID: 8291699 DOI: 10.1097/00000542-199401000-00018
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = 3
        self.age_upper = 11
        self.weight_lower = 15
        self.weight_upper = 61
        self.pmid = "8291699"
        self.doi = "10.1097/00000542-199401000-00018"
        self.validate_anthropometric_values()

        self.v1 = weight * 0.41
        self.v2 = weight * 0.78 + 3.1 * age - 16
        self.v3 = weight * 6.9
        self.cl1 = weight * 0.035
        self.cl2 = weight * 0.077
        self.cl3 = weight * 0.026

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k13 = self.cl3 / self.v1
        self.k21 = self.cl2 / self.v2
        self.k31 = self.cl3 / self.v3
        self.ke0 = 0.41


class Eleveld(Model):
    """Eleveld class holds pharmacokinetic parameters for the Eleveld propofol
    model.

    Reference: PMID: 29661412 DOI: 10.1016/j.bja.2018.01.018
    """

    def __init__(self, sex: int, age: float, weight: float, height: float,
                 opiates_coadministered: bool):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = -1
        self.age_upper = -1
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "29661412"
        self.doi = "10.1016/j.bja.2018.01.018"
        self.validate_anthropometric_values()

        # Reference
        sex_ref = 0
        age_ref = 35
        weight_ref = 70
        height_ref = 170

        # Theta constants
        theta_1 = 6.28  # V1ref L
        theta_2 = 25.5  # V2ref L
        theta_3 = 273  # V3ref L
        theta_4 = 1.79  # CLref (male) L/min
        theta_5 = 1.83  # Q2ref L/min
        theta_6 = 1.11  # Q3ref L/min
        theta_7 = 0.191  # Typical residual error
        theta_8 = 42.3  # CL maturation E50 weeks
        theta_9 = 9.06  # CL maturation slope
        theta_10 = -0.0156  # Smaller V2 with age
        theta_11 = -0.00286  # Lower CL with age
        theta_12 = 33.6  # Weight for 50% of maximal V1 Kg
        theta_13 = -0.0138  # Smaller V3 with age
        theta_14 = 68.3  # Maturation of Q3 weeks
        theta_15 = 2.1  # CLref (female) L$min1
        theta_16 = 1.3  # Higher Q2 for maturation of Q3
        theta_17 = 1.42  # V1 venous samples (children)
        theta_18 = 0.68  # Higher Q2 venous samples

        # Post menstrual age
        pma = age * 52 + 40
        pma_ref = age_ref * 52 + 40

        def ageing(x, age):
            return exp(x * (age - age_ref))

        def sigmoid(x, e50, y):
            return (x ** y) / ((x ** y) + (e50 ** y))

        def central(x):
            return sigmoid(x, theta_12, 1)

        def opiates(x, present):
            if present:
                return exp(x * age)
            else:
                return 1

        # cl1 maturation
        cl1_mat = sigmoid(pma, theta_8, theta_9)
        cl1_mat_ref = sigmoid(pma_ref, theta_8, theta_9)

        # cl3 maturation
        cl3_mat = sigmoid(pma, theta_14, 1)
        cl3_mat_ref = sigmoid(pma_ref, theta_14, 1)

        # fat free mass
        ffm = ffm_alsallami(sex, age, weight, height)
        ffm_ref = ffm_alsallami(sex_ref, age_ref, weight_ref, height_ref)

        self.v1 = theta_1 * (central(weight) / central(weight_ref))
        self.v2 = theta_2 * (weight / weight_ref) * ageing(theta_10, age)
        self.v3 = theta_3 * (ffm / ffm_ref) * opiates(theta_13,
                                                      opiates_coadministered)

        if sex == 0:
            self.cl1 = theta_4 * ((weight / weight_ref) ** 0.75) * \
                       (cl1_mat / cl1_mat_ref) * \
                       opiates(theta_11, opiates_coadministered)
        elif sex == 1:
            self.cl1 = theta_15 * (weight / weight_ref) ** 0.75 * \
                       (cl1_mat / cl1_mat_ref) * \
                       opiates(theta_11, opiates_coadministered)

        self.cl2 = theta_5 * (self.v2 / theta_2) ** 0.75 * \
                    (1 + theta_16 * (1 - cl3_mat))

        self.cl3 = theta_6 * (self.v3 / theta_3) ** 0.75 * \
                   (cl3_mat / cl3_mat_ref)

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k13 = self.cl3 / self.v1
        self.k21 = self.cl2 / self.v2
        self.k31 = self.cl3 / self.v3
        self.ke0 = 0.146 * ((weight / weight_ref) ** -0.25)

        self.ce50 = 3.08 * ageing(-0.00635, age)


class Short(Model):
    """Short class holds pharmacokinetic parameters for the Short
    propofol model; uses Keo from Eleveld

    Reference: PMID: 8130049 DOI: 10.1093/bja/72.3.302
    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = 4
        self.age_upper = 7
        self.weight_lower = 15
        self.weight_upper = 22
        self.pmid = "8130049"
        self.doi = "10.1093/bja/72.3.302"
        self.validate_anthropometric_values()

        self.v1 = 0.432 * weight

        self.k10 = 0.0967
        self.k12 = 0.1413
        self.k13 = 0.1092
        self.k21 = 0.0392
        self.k31 = 0.0049
        self.ke0 = 0.146 * ((weight / 70) ** -0.25)


class Schuttler(Model):
    """Schuttler class holds pharmacokinetic parameters for the Schuttler
       propofol model; uses Keo from Eleveld

    Reference: PMID:  DOI: 10.1097/00000542-200003000-00017
    """

    def __init__(self, sex: int, age: float, weight: float, height: float,
                 bolus_data: bool = False, venous_data: bool = False):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = 2
        self.age_upper = 88
        self.weight_lower = 2
        self.weight_upper = 88
        self.pmid = "10719952"
        self.doi = "10.1097/00000542-200003000-00017"
        self.validate_anthropometric_values()

        # Theta constants
        theta_1 = 1.44
        theta_2 = 9.3
        theta_3 = 2.25
        theta_4 = 44.2
        theta_5 = 0.92
        theta_6 = 266
        theta_7 = 0.75
        theta_8 = 0.62
        theta_9 = 0.61
        theta_10 = 0.045
        theta_11 = 0.55
        theta_12 = 0.71
        theta_13 = 20.39
        theta_14 = 20.40
        theta_15 = 1.61
        theta_16 = 2.02
        theta_17 = 0.73
        theta_18 = 20.48

        if bolus_data:
            bol = 1
        else:
            bol = 0

        if venous_data:
            ven = 1
        else:
            ven = 0

        self.v1 = theta_2 * (weight / 70) ** theta_12 * \
                  (age / 30) ** theta_13 * (1 + bol * theta_15)
        self.v2 = theta_4 * (weight / 70) ** theta_9 * \
                  (1 + bol * theta_17)
        self.v3 = theta_6

        if age <= 60:
            self.cl1 = theta_1 * (weight / 70) ** theta_7
        else:
            self.cl1 = theta_1 * (weight / 70) ** theta_7 - \
                       (age - 60) * theta_10
        self.cl2 = theta_3 * (weight / 70) ** theta_8 * \
                   (1 + ven * theta_14) * (1 + bol * theta_16)
        self.cl3 = theta_5 * (weight / 70) ** theta_11 * \
                   (1 + bol * theta_18)

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k13 = self.cl3 / self.v1
        self.k21 = self.cl2 / self.v2
        self.k31 = self.cl3 / self.v3
        self.ke0 = 0.146 * ((weight / 70) ** -0.25)
