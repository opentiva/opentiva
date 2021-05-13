"""
opentiva.helpers
================

This module contains the helper functions for model calculations.
"""


def lbm_dubois(sex: int, weight: float, height: float) -> float:
    """Returns lean body mass using DuBois method

    Parameters
    ----------
    sex
        0 for male or 1 for female
    weight
        weight in kg
    height
        height in cm

    Returns
    -------
    float
        lean body mass
    """
    if sex == 0:
        return 1.1 * weight - 128 * (weight / height) ** 2
    elif sex == 1:
        return 1.07 * weight - 148 * (weight / height) ** 2


def body_mass_index(weight: float, height: float) -> float:
    """Returns body mass index

    Parameters
    ----------
    weight
        weight in kg
    height
        height in cm

    Returns
    -------
    float
        body mass index
    """
    return round(weight / (height / 100) ** 2, 1)


def ffm_janmahasation(sex: int, weight: float, height: float) -> float:
    """ Returns fat free mass using Janmahasation method

    Parameters
    ----------
    sex
        0 for male or 1 for female
    weight
        weight in kg
    height
        height in cm

    Returns
    -------
    float
        fat free mass
    """

    bmi = body_mass_index(weight, height)

    if sex == 0:
        ffm = (9270 * weight) / (6680 + 216 * bmi)
    elif sex == 1:
        ffm = (9270 * weight) / (8780 + 244 * bmi)

    return ffm


def ffm_alsallami(sex: int, age: float, weight: float, height: float) -> float:
    """ Method returns fat free mass using Alsallami method

    Parameters
    ----------
    sex
        0 for male or 1 for female
    age
        age in years
    weight
        weight in kg
    height
        height in cm

    Returns
    -------
    float
        fat free mass
    """

    if sex == 0:
        ffm = (0.88 + (
            (1 - 0.8) / (1 + (age / 13.4) ** -12.7)
        )) * ffm_janmahasation(sex, weight, height)
    elif sex == 1:
        ffm = (1.11 + (
            (1 - 1.11) / (1 + (age / 7.1) ** -1.1)
        )) * ffm_janmahasation(sex, weight, height)

    return ffm


def crcl_cockcroft_gault(sex: int, age: float, weight: float, height: float,
                         creatinine: float) -> float:
    """ Method returns creatinine clearance using the Cockcroft-Gault  method

    Parameters
    ----------
    sex
        0 for male or 1 for female
    age
        age in years
    weight
        weight in kg
    height
        height in cm
    creatinine
        serum creatinine value in umol/L

    Returns
    -------
    float
        creatinine clearance
    """

    if sex == 0:
        ibw_f = 50
        sex_f = 1
    else:
        ibw_f = 45.5
        sex_f = 0.85

    if height - 152.4 > 0:
        height_f = height - 152.4
    else:
        height_f = 0

    ibw = ibw_f + (2.3 * height_f)

    bmi = body_mass_index(weight, height)
    if bmi < 18.5:
        w = weight
    elif bmi >= 18.5 and bmi < 25:
        w = ibw
    elif bmi >= 25:
        w = ibw + 0.4 * (weight - ibw)

    crcl = ((140 - age) * w * sex_f) / (creatinine * (72 / 88.42))

    return crcl


def crcl_schwartz(height: float, creatinine: float) -> float:
    """ Method returns creatinine clearance using the Schwartz method

    Parameters
    ----------
    height
        height in cm
    creatinine
        serum creatinine value in umol/L

    Returns
    -------
    float
        creatinine clearance
    """

    k = 0.413
    crcl = k * height / (creatinine / 88.42)

    return crcl


def bsa_dubois(weight: float, height: float) -> float:
    """ Method returns body surface area using Debois method

    Parameters
    ----------
    weight
        weight in kg
    height
        height in cm

    Returns
    -------
    float
        body surface area
    """

    return 0.007184 * height ** 0.725 * weight ** 0.425
