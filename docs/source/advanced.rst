Advanced
========

User-defined infusions
----------------------

User defined infusion are where the start time, dose and duration of an 
infusion are inputted.
These can be used with or without target controlled infusions.

.. warning::
   Using an user defined infusion alongside target controlled infusions may result in overshoot of the target

Example
~~~~~~~

.. code:: python
  
   import opentiva
   
   p1 = opentiva.pump.Pump(...)

   p1.add_infusion(start=300, dose=0.042, duration=(60*60))

The above example adds an infusion starting at 5 minutes (300 seconds) at a dose of
0.042 mg / sec over 1 hour (3600 seconds) to give a total infused dose of
151 mg.

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description
   * - start
     - start time of infusion in seconds
   * - dose
     - dose per second (e.g. mg / second) of infusion; unit as specified in model.concentration_unit
   * - duration
     - total time of infusion in seconds

.. note::
   Bolus doses can be considered as an infusion with a small duration 
   e.g. dose 150 over duration 1 second 


Decrement time
--------------

Decrement time is the time in seconds to reach a target from a start time
if all infusions were stopped.

Example
~~~~~~~

.. code:: python
  
   import opentiva
   
   p1 = opentiva.pump.Pump(...)
   p1.add_target(start=0, target=4, duration=10, effect=False)

   decrement_time = p1.decrement_time(start=300, effect=False,
                                      target=1)

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description
   * - start
     - start time to calculate decrement time from in seconds
   * - effect
     - true for effect site decrement time or false for plasma site
       decrement time
   * - target
     - target concentration of plasma or effect concentration to
       decrement to

.. note::
   Due to the exponential nature of the decrement time a concentration of 0 
   occurs at infinity. To avoid infinite loops a selected decrement target of 
   0 will approximate to 0.1.

Exporting data
--------------

Data can be exported by directly calling the Pump class object or by using a helper function:

Direct calling
~~~~~~~~~~~~~~

**target_concentrations**

.. code:: python
  
   import opentiva
   
   p1 = opentiva.pump.Pump(...)

   output = p1.target_concentrations

2d array of targets with each row containing:

* start - start time of target in seconds
* target - target of plasma or effect concentration
* duration - time to achieve target level in seconds
* effect - true for effect site targetting or false for plasma site targetting
* cp_limit - multiplied to the target to get maximum increase in plasma
  concentration
* cp_limit_duration - time to achieve cp_limit in seconds
* ce_bolus_only - true will target the effect using the 
  'original' method / bolus only method

.. note::
   In original effect site targetting after running generate_infusions() or 
   run() functions the cp_limit for that target will be replaced with the 
   calculated cp_limit value.

**infusion_list**

.. code:: python

   p1.generate_infusions()
   output = p1.infusion_list

2d array of infusions calculated from provided targets and including any user 
defined infusions. With each row containing:

* start time of infusion in seconds
* dose of infusion over 1 second
* duration of infusion in seconds
* end time of infusion in seconds

**user_infusion_list**

.. code:: python
  
   output = p1.user_infusion_list

2d array of infusions entered by the user with each row containing:

* start time of infusion in seconds
* dose of infusion over 1 second
* duration of infusion in seconds
* end time of infusion in seconds

Helper functions
~~~~~~~~~~~~~~~~

**generate_rates_array**

.. code:: python
  
   output = p1.generate_rates_array()

Turns the infusions from the infusion_list array into a ml/hr array; 
with each row containing:

* time of rate change (seconds)
* rate (ml/hr)

**generate_dose_weight_array**

.. code:: python

   output_min = p1.generate_dose_weight_array('min')
   output_hr = p1.generate_dose_weight_array('hr')

Turns the infusions from the infusion_list array into dose/weight 
(if the infusion time is below the pump's bolus time) or dose/weight/time
array if not; with each row containing:

* time of rate change (seconds);
* bolus dose expressed as dose/weight (e.g. mg/kg) or
  dose/weight/time (e.g. mg/kg/hr);
* true if bolus dose/weight (e.g. mg/kg) or
  false if infusion dose/weight/time (e.g. mg/kg/hr)

The time interval for dose/weight/time unit can be either:

* min - str - for dose/weight/minute (default)
* hr - str - for dose/weight/hour

**generate_targets_array**

.. code:: python

   output = p1.generate_targets_array()

Turns the target_concentrations array into a simple targets array 
with each row containing:

* time of target change (seconds)
* target concentration

Add new drug models
-------------------

New drug models can be easily added but need to follow a default format. 
They should be a python class with the parameters and
attributes shown in the minimal example. The full example is the recommeded
format over the minimal. 

For validation of the entered sex, age, weight and height the model should 
inherit the parent class from opentiva.model.

The class can then be imported to the pump. 

Take a new module `newdrug` with model as class `Model`:

.. code:: python

   import opentiva
   from . import newdrug

   drug_model = newdrug.Model(sex=sex, age=age,
                              weight=weight, height=height)

   p1 = opentiva.pump.Pump(model=drug_model, drug_concentration=10,
                           end_time=(60*60))

Minimal example
~~~~~~~~~~~~~~~

`./drug.py`

.. code:: python
   
   class MinimalModel():

       def __init__(self, sex: int, age: float, weight: float, height: float):
           self.weight = weight
           self.compartments = 3

           self.v1 = 14.3  # unit litres

           self.k10 = 0.0645  # unit /min
           self.k12 = 0.1086  # unit /min
           self.k13 = 0.0229  # unit /min
           self.k21 = 0.0245  # unit /min
           self.k31 = 0.0014  # unit /min
           self.ke0 = 0.17559  # unit /min

Full example
~~~~~~~~~~~~

Annotated example of the propofol marsh diprifusor model.

.. code:: python

   import opentiva.model as Model
     # parent class to validate class instance variables
   import opentiva.validation as validation  
     # functions to validate anthopometric values


      class MarshDiprifusor(Model):
          """MarshDiprifusor class holds pharmacokinetic parameters for the
             Diprifusor Marsh propofol model with Keo 0.26.
             Reference: PMID: 1859758 DOI: 10.1093/bja/67.1.41
          """
   
          def __init__(self, sex: int, age: float, weight: float, height: float):
              super().__init__(sex, age, weight, height)

              # sex : 0 for male or 1 for female
              # age : in years
              # weight : in kg
              # height : in cm

              self.compartments = 3
                # number of compartments to model; either 2 or 3
              self.concentration_unit = "mg/ml"
                # drug concentration unit
              self.target_unit = "ug/ml"
                # target concentration unit
              self.age_lower = 16
                # lower age limit of model; -1 if no limit
              self.age_upper = -1
                # upper age limit of model; -1 if no limit
              self.weight_lower = -1
                # lower weight limit of model; -1 if no limit
              self.weight_upper = 150
                # upper weight limit of model; -1 if no limit
              self.pmid = "1859758"
                # pubmed ID of model's reference
              self.doi = "10.1093/bja/67.1.41"
                # digital Object Identifier (DOI) of model's reference
              self.validate_anthropometric_values()
                # Funtion imported from the parent Model class and validates 
                # the input age and weight is between the specified values of 
                # the model. If not a warning is appended to the self.warning 
                # string and to sys.stderr.
   
              self.v1 = 0.228 * weight
                # volume of central compartment in litres
              self.v2 = 0.463 * weight
                # volume of fast compartment in litres
              self.v3 = 2.893 * weight
                # volume of slow compartment in litres
   
      #       self.cl1 =
                # clearence of compartment 1 in litres/min
      #       self.cl2 =
                # clearence of compartment 2 in litres/min
      #       self.cl3 =
                # clearence of compartment 3 in litres/min

              self.k10 = 0.119
                # equilibrium rate constant from compartment 1 to 0 /min
                # also = self.cl1 / self.v1
              self.k12 = 0.112
                # equilibrium rate constant from compartment 1 to 2 /min
                # also = self.cl2 / self.v1
              self.k13 = 0.0419
                # equilibrium rate constant from compartment 1 to 3 /min
                # also = self.cl3 / self.v1
              self.k21 = 0.055
                # equilibrium rate constant from compartment 2 to 1 /min
                # also = self.cl2 / self.v2
              self.k31 = 0.0033
                # equilibrium rate constant from compartment 3 to 1 /min
                # also = self.cl3 / self.v3
              self.ke0 = 0.26
                # effect compartment equilibrium rate constant /min
   
      #       self.k20 = 
                # equilibrium rate constant from compartment 2 to 0 /min
                # optional parameter for two compartment modelling


Deriving Ke0
------------

Ke0 can be derived from a PKPD model using the :math:`t_{peak}` method. 
Using a preexisiting model or :ref:`a new model<Add new drug models>` 
(if using a new model and ke0 is unknown then self.ke0 must still be 
declared within the module e.g. self.ke0 = 0)

Details of the theory surrounding this can be found :ref:`here<Ke0 'tpeak' method>`.

.. code:: python

    import opentiva

    sex = 0  # male
    age = 30  # years
    weight = 70  # kg
    height = 170  # cm

    propofol_marsh = opentiva.propofol.MarshDiprifusor(sex=sex, age=age,
                                                       weight=weight, 
                                                       height=height)
    
    dose = 1 # bolus dose in mg
    tpeak = 236 # time to peak effect site concentration in seconds
    ce_tpeak = 0.25831 # effect site concentration at tpeak

    ke0 = propofol_marsh.ke0_tpeak_method(dose=dose, tpeak=tpeak, 
                                          ce_tpeak=ce_tpeak)

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Parameter
     - Description
   * - dose
     - total bolus dose of drug given at time 0 over 1 second
   * - tpeak
     - time in seconds of peak effect
   * - ce_tpeak
     - effect site concentration at tpeak

.. note::
   ke0_tpeak_method returns the ke0 as /seconds.
