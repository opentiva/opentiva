from numpy import exp

from . helpers import bsa_dubois, crcl_cockcroft_gault, crcl_schwartz
from . model import Model

"""
opentiva.rocuronium
===================

This module contains the classes for rocuronium models.
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
asian
   bool, true if race is asian, Kleijn only
creatinine
   creatinine value, Kleijn only
sevoflurane
   bool, true if sevoflurane and false if TIVA, Kleijn only

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


class Kleijn(Model):
    """Kleijn class holds pharmacokinetic parameters for the Kleijn
    rocuronium model with no sevoflurane.

    Reference: PMID: 21535448 DOI: 10.1111/j.1365-2125.2011.04000.x

    """

    def __init__(self, sex: int, age: float, weight: float, height: float,
                 creatinine: float, sevoflurane: bool, asian: bool = False):
        super().__init__(sex, age, weight, height)

        self.compartments = 2
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = -1
        self.age_upper = -1
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "21535448"
        self.doi = "10.1111/j.1365-2125.2011.04000.x"
        self.validate_anthropometric_values()

        if age >= 18:
            crcl = crcl_cockcroft_gault(sex, age, weight, height, creatinine)
        else:
            crcl = crcl_schwartz(height, creatinine)
            crcl *= bsa_dubois(weight, height)  # denomalize using BSA

        v1_cr = exp(-0.00143 * (crcl - 119))
        self.v1 = v1_cr * 4.73 * (weight / 70)

        cl_age = 1 + -0.00678 * (age - 43)
        self.cl1 = cl_age * 0.269 * (weight / 70) ** 0.75

        v2_age = exp(0.00613 * (age - 43))
        self.v2 = v2_age * 6.76 * (weight / 70)

        if asian:
            q2_rac = 1 + -0.212
        else:
            q2_rac = 1

        self.cl2 = q2_rac * 0.279 * (weight / 70) ** 0.75

        self.k10 = self.cl1 / self.v1
        self.k12 = self.cl2 / self.v1
        self.k21 = self.cl2 / self.v2

        if sevoflurane:
            ke0_sev = 1 + -0.567
        else:
            ke0_sev = 1

        self.ke0 = ke0_sev * 0.134 * (weight / 70) ** -0.25

    @property
    def creatinine(self):
        return self._creatinine

    @creatinine.setter
    def creatinine(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError("Expected: float or int")
        if value <= 0:
            raise ValueError("Must be greater than 0")
        self._creatinine = value


class Wierda(Model):
    """Wierda class holds pharmacokinetic parameters for the Wierda
    rocuronium model.

    Reference: PMID: 1829656 DOI: 10.1007/BF03007578

    Ke0: PMID: 30099599 DOI: 10.1007/s00540-018-2543-3

    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 3
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = 18
        self.age_upper = 60
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "1829656"
        self.doi = "10.1007/BF03007578"
        self.validate_anthropometric_values()

        self.v1 = 0.044 * weight

        self.k10 = 0.1
        self.k12 = 0.21
        self.k13 = 0.028
        self.k21 = 0.13
        self.k31 = 0.01

        self.ke0 = 0.1 + (age - 50) * -0.000725


class Woloszczuk(Model):
    """Woloszczuk class holds pharmacokinetic parameters for the Woloszczuk
    rocuronium model.

    Reference: PMID: 16879519 DOI: 10.1111/j.1460-9592.2005.01840.x

    """

    def __init__(self, sex: int, age: float, weight: float, height: float):
        super().__init__(sex, age, weight, height)

        self.compartments = 2
        self.concentration_unit = "mg/ml"
        self.target_unit = "ug/ml"
        self.age_lower = 3
        self.age_upper = 11
        self.weight_lower = -1
        self.weight_upper = -1
        self.pmid = "16879519"
        self.doi = "10.1111/j.1460-9592.2005.01840.x"
        self.validate_anthropometric_values()

        self.v1 = 28.36 * weight / 1000

        self.k10 = 0.03
        self.k12 = 0.125
        self.k21 = 0.003
        self.ke0 = 0.1
