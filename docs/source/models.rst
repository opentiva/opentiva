Drug Models
===========

Attributes 
----------

The drug models all share the same availiable attributes shown below. These
can be returned by interacting with the model object.

.. code:: python

    import opentiva

    sex = 0  # male
    age = 30  # years
    weight = 70  # kg
    height = 170  # cm

    remifentanil_minto = opentiva.remifentanil.Minto(sex=sex, age=age,
                                                     weight=weight, height=height)

    print(remifentanil_minto.compartments)
        # Number of compartments that the model uses
    print(remifentanil_minto.k12)
        # Equilibrium rate constant from compartment 1 to 2
    print(remifentanil_minto.ke0)
        # Effect compartment equilibrium rate constant

.. list-table::
   :widths: 33 33 33
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - compartments
     - int
     - number of compartments to model; 1, 2 or 3
   * - concentration_unit
     - str
     - drug concentration unit
   * - target_unit 
     - str
     - target concentration unit
   * - age_lower
     - int
     - lower age limit of model; -1 if no limit
   * - age_upper
     - int
     - upper age limit of model; -1 if no limit
   * - weight_lower
     - int
     - lower weight limit of model; -1 if no limit
   * - weight_upper
     - int
     - upper weight limit of model; -1 if no limit
   * - pmid
     - str
     - Pubmed ID of model's reference
   * - doi
     - str
     - Digital Object Identifier (DOI) of model's reference
   * - warning
     - str
     - Warnings relating to non-validated anthropometric values
   * - v1
     - float
     - volume of central compartment (litres)
   * - k10
     - float
     - equilibrium rate constant from compartment 1 to 0 (/min)
   * - k12
     - float
     - equilibrium rate constant from compartment 1 to 2 (/min)
   * - k13
     - float
     - equilibrium rate constant from compartment 1 to 3 (/min)
   * - k21
     - float
     - equilibrium rate constant from compartment 2 to 1 (/min)
   * - k31
     - float
     - equilibrium rate constant from compartment 3 to 1 (/min)
   * - ke0
     - float
     - effect compartment equilibrium rate constant (/min)


Hypnotics
---------

* Dexmedetomidine
   * Dyck
   * Hannivoort
   * PerezGuille
   * Rolle
* Ketamine
   * Clements250
   * Domino
   * Herd
   * Hijazi
   * Hornik
   * Klamp
* Midazolam
   * Albrecht
* Propofol
   * Eleveld (no opiates)
   * Eleveld (opiates)
   * Kataria
   * Marsh Diprifusor
   * Marsh Modified
   * Paedfusor
   * Schnider
   * Schuttler
   * Short
* Thiopental
   *  Stanski

Opioids
-------

* Alfentanil
   * Maitre
   * Goresky
   * Scott
* Fentanyl
   * Ginsberg
   * Scott
   * Shafer
   * ShaferW80
* Morphine
   * Sarton
* Remifentanil
   * Eleveld
   * Kim
   * Minto
   * Rigby-Jones
* Sufentanil
   * Gepts
   * Greely

Neuromuscular blockers
----------------------

* Atracurium
   *  Fisher
   *  Marathe
* Cisatracurium
   *  Bergeron
   *  Imbeault
   *  Tran
* Rocuronium
   *  Kleijn Sevo
   *  Kleijn TIVA
   *  Wierda
   *  Woloszczuk
* Vecuroniumg
   *  Caldwell
   *  Wierda

