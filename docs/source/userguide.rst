Userguide
=========

Overview
--------

There are two main steps in configuring a simulation.

#. :ref:`Drug and model selection<Drug and model selection>`
#. :ref:`Pump and target configuration<Pump configuration>`


Drug and model selection
------------------------

All the availiable drug models follow the same import format. 

.. code:: python

    import opentiva

    drug_model = opentiva.drug.Model(sex, age, weight, height)

The parameters are:

.. list-table:: 
   :widths: 33 33 33
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - Sex
     - int
     - 0 for male or 1 for female
   * - Age
     - int
     - in years
   * - Weight
     - float
     - in kg
   * - Height
     - float
     - in cm

List of availiable models can be found :ref:`here<Drug models>`.

.. note::
   Some models require extra parameters; for example, propofol Eleveld's 
   opiates_coadministered bool (True if opiate coadministration) or rocuronium 
   Kleijn's sevoflurane, asian and creatinine parameters.

Examples
~~~~~~~~

**Importing:**

.. code:: python

    import opentiva

    sex = 0  # male
    age = 30  # years
    weight = 70  # kg
    height = 170  # cm

    propofol_schnider = opentiva.propofol.Schnider(sex=sex, age=age,
                                                   weight=weight, height=height)
    propofol_eleveld = opentiva.propofol.EleveldWithoutOpiates(sex=sex, age=age,
                                                               weight=weight, height=height)
    remifentanil_minto = opentiva.remifentanil.Minto(sex=sex, age=age, 
                                                     weight=weight, height=height)

**Returning a model's pharmokinetic/ pharmacodynamic values:**

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

All the availiable exportable pharmokinetic/ pharmacodynamic values can be 
found :ref:`here<Attributes>`.

**Returning a model's warnings:**

Warnings will be generated when a models validated anthropometric values are 
not met. These will written to sys.stderr and returned as a string.

.. code:: python

    import opentiva

    sex = 1  # female
    age = 5  # years
    weight = 18  # kg
    height = 109  # cm

    propofol_marsh = opentiva.propofol.MarshDiprifusor(sex=sex, age=age,
                                                       weight=weight, height=height)
    warning = propofol_marsh.warning
    print(warning)

Pump configuration
------------------

The pump class has a range of required and optional settings that determine 
how the simulation runs.

A minimal example of the required parameters is:

.. code:: python

    import opentiva
    
    propofol_eleveld = opentiva.propofol.Eleveld(sex=0, age=30, 
                                                 weight=70, height=170,
                                                 opiates_coadministered=False)
    
    p1 = opentiva.pump.Pump(model=propofol_eleveld, drug_concentration=10,
                            end_time=(60*60))

.. list-table:: 
   :widths: 33 33 33
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - model
     - Class object
     - opentiva pharmacokinetic/ pharmacodynamic model
   * - drug_concentration
     - float
     - concentration of infusion drug (using units outlined in model; `model.concentration_unit`)
   * - end_time
     - int
     - duration of simulation in seconds


The optional settings relate to defaults for effect site targetting, maintenance
infusions and characteristics of the simulated pump. 

Full API reference can be found :ref:`here<modindex>`.

The effect site targetting parameters (cp_limit and cp_limit_duration) are 
covered :ref:`here<Effect site targetting>`.

Maintenance infusions
~~~~~~~~~~~~~~~~~~~~~

The maintenance infusion is the infusion required to offset the clearance and 
distribution losses to maintain a steady state plasma and effect site 
concentration.

As the time in a steady state increases the exponential nature of the two/ 
three compartment model means the changes between the amount of drug being 
infused by the each subsequent maintenance infusions decreases.

**maintenance_infusion_duration** int *default 300 seconds*

Time in seconds of the duration of each maintenance infusion

**maintenance_infusion_multiplier** float *default 2*

The duration of each subsequent maintenance infusion is multiplied by this, 
extending the time between each maintenance infusion calculation.

This reduces the number of infusions required and computing load, but may 
result in periods of time below a target.


Pump characteristics
~~~~~~~~~~~~~~~~~~~~

**max_infusion_rate** float *default 1200 ml/hr*

Maximum rate in ml/hr that the pump can deliver an infusion. 
If a calculated infusion rate is higher than this it will be limited to this rate. 
This can increase the time to achieve a target plasma and effect site concentration.

**bolus_time** int *default 20 seconds*

If the duration of an infusion is less than this time the dose will be 
considered as being given as a bolus and the rate of infusion will not be 
limited by the set maximum infusion rate.

Full example
~~~~~~~~~~~~

.. code:: python

    import opentiva
    
    propofol_paedfusor = opentiva.propofol.Paedfusor(sex=0, age=5, 
                                                     weight=18, height=109)
    
    p1 = pump.Pump(model=propofol_paedfusor,
                   drug_concentration=10,  # propofol concentration of 10 mg/ml
                   end_time=(2*60*60),  # total duration of simulation is 2 hours
                   maintenance_infusion_duration=10,  # first maintenance infusion is 10 seconds in duration
                   maintenance_infusion_multiplier=2,  # the duration of each subsequent maintenance infusion is doubled
                   cp_limit = 1.2,  # maximum plasma concentration during effect targetting is 1.2 x the target
                   cp_limit_duration = 10,  # duration to achieve cp_limit is 10 seconds
                   max_infusion_rate = 1200,  # maximum rate of pump is 1200 ml/hr
                   bolus_time = 20)  # any infusion < 20 seconds in duration is considered a bolus dose


Target configuration
--------------------

Once the drug/ model and pump have been initialized target concentrations can be added. 
These can target either the plasma or effect site.

This is achieved by using the `add_target` function.

.. code:: python

    import opentiva

    p1 = opentiva.pump.Pump(...)

    p1.add_target(...)

Plasma site targetting
~~~~~~~~~~~~~~~~~~~~~~

In plasma site targetting opentiva calculates the infusion required to reach 
the plasma target concentration after a specified duration.

.. code:: python

    p1.add_target(start=0, target=4, duration=10, effect=False)

The above example will reach a plasma concentration of 4 after 10 seconds 
starting from 0 seconds.

.. list-table:: 
   :widths: 33 33 33
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - start
     - int
     - start time of target in seconds
   * - target
     - float
     - target of plasma or effect concentration
   * - duration
     - int
     - time to achieve target level in seconds
   * - effect
     - bool: default *True*
     - true for effect site targetting or false for plasma site targetting



Effect site targetting
~~~~~~~~~~~~~~~~~~~~~~

In effect site targetting there are two methods availiable to calculate 
the infusions required to reach the effect site concentration.

**Original**

Original method described by 
`Shafer et al <(https://link.springer.com/article/10.1007/BF01070999)>`_.

An infusion is given over a duration *(Cp_limit_duration)* which causes 
an overshoot in plasma concentration. 
This overshoot decreases the time to reach a specified 
effect site target *(Cetarget)* compared to plasma targetting mode.

.. image:: ../images/original.png
  :width: 600
  :alt: original effect site targetting

Example:

.. code:: python

    p1.add_target(start=0, target=4, duration=10, effect=True,
                  cp_limit_duration=20, ce_bolus_only=True)

The above example will calculate the bolus dose over 20 seconds to achieve 
a rapid rise to an effect site of 4 with minimal effect site overshoot.

.. list-table:: 
   :widths: 33 33 33
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - cp_limit_duration
     - int
     - duration in seconds to achieve the maximal plasma concentration
   * - ce_bolus_only
     - bool
     - true for Original method (bolus only) effect targetting or false for revised method

**Revised**

Revised effect targetting aims to decrease the overshoot in plasma 
concentration; 
based on `Van Poucke et al <https://ieeexplore.ieee.org/document/1344189>`_. 

To achieve a target effect site *(Cetarget)* an initial infusion to reach a 
maximal plasma target *(Cpmax)* is given. 
This plasma target is set by the multiplying the target effect 
concentration by a limit value *(Cplimit)*. 
Once at this maximum plasma concentration a maintenance infusion is started 
to maintain steady state over time *(Tinf)*. 
At a certain point in time stopping the maintenance infusion will result in 
the plasma concentration decrementing over time *(Tcoast)* to the effect target 
concentration which will be met at the same time as the up trending effect 
site concentration.

.. image:: ../images/revised.png
  :width: 600
  :alt: revised effect site targetting

Example:

.. code:: python

   p1.add_target(start=0, target=4, duration=10, effect=True,
                 cp_limit=1.5, cp_limit_duration=20,
                 ce_bolus_only=True)

The above example will increase to a plasma concentration of 6 over 20 seconds.
This is then maintained at 6 over a calculated duration then stopped. After stopping 
the decrement of the plasma concentration 
meets the increasing effect site concentration at the target of 4.

.. list-table:: 
   :widths: 33 33 33
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - cp_limit
     - float
     - multiplied to the target to get maximum increase in plasma concentration

.. note::
   **Optimal value for cp_limit:** The maximum value of cp_limit can be found by
   first running the simulation with the target reached by the original method 
   (ce_bolus_only=True). Above this calculated cp_limit an overshoot in effect
   site concentration would occur. After generating infusions (generate_infusions())
   or running the simulation (run()) the calculated cp_limit can be obtained from the
   :ref:`target_concentration<Direct calling>` array. Then the cp_limit for the 
   revised targetting can be set a value less than this.

**Duration in effect targetting**

The *duration* parameter will only have an influence on effect site 
targetting if it is longer then the minimum time it takes to reach a target.

For example, if it takes 250 seconds for a effect site target to be reached 
by the original method.

- A *duration* of 10 seconds will be ignored.
- A *duration* of 300 seconds will change the infusions so that the effect site 
  target is reached at approximately 300 seconds.

Maintenance infusions
~~~~~~~~~~~~~~~~~~~~~

By default opentiva will calculate the maintenance infusions required to 
maintain steady state from when the target is reached to the next set target. 

To disable this set *maintenance_infusions* as False:

.. code:: python

   p1.add_target(start=0, target=4, duration=10, effect=True,
                 cp_limit=1.5, cp_limit_duration=20, ce_bolus_only=True,
                 maintenance_infusions=False)

The above example does the same as the revised effect site targetting example 
but will not maintain a steady state after the target is reached.

.. list-table:: 
   :widths: 33 33 33
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - maintenance_infusions
     - bool: default *True*
     - true will calculate the infusions between the time the target is 
       reached and the next target to maintain steady state

Run simulation
--------------

The simulation can then be run by using the run function:

.. code:: python

   output = p1.run()

This returns a 2d numpy array of:

.. list-table:: 
   :widths: 50 50
   :header-rows: 1

   * - Column
     - Description
   * - 1
     - time in seconds of concentration
   * - 2
     - plasma concentration   
   * - 3
     - effect site concentration

To generate the infusions required to meet the targets and not calculate 
the plasma and effect site concentrations the following function can be used.

.. code:: python

   p1.generate_infusions()

See :ref:`Exporting data` on how to export the infusions.

Putting it all together
~~~~~~~~~~~~~~~~~~~~~~~

A full example is show below:

.. code:: python

    import opentiva

    adult35 = opentiva.propofol.Eleveld(sex=0, age=35, 
                                        weight=70, height=170, 
                                        opiates_coadministered=False)

    p1 = opentiva.pump.Pump(model=adult35,
                           drug_concentration=10,
                           end_time=(1*60*60),
                           maintenance_infusion_duration=10,
                           maintenance_infusion_multiplier=2,
                           cp_limit = 1.2,
                           cp_limit_duration = 10,
                           max_infusion_rate = 1200,
                           bolus_time = 20)

    # Target Ce 3 using revised method at time 0
    p1.add_target(start=0,
                  target=3,
                  duration=10,
                  effect=True,
                  cp_limit=5.0,
                  cp_limit_duration=10,
                  ce_bolus_only=False)

    # Reduce to Ce 2.5 at 5 minutes (300 seconds)
    p1.add_target(start=300,
                  target=2.5,
                  duration=10,
                  effect=True)

    # Target Cp 4 at time 10 minutes (600 seconds) over 1 minute
    p1.add_target(start=1200,
                  target=4,
                  duration=60,
                  effect=False)

    # Reduce to Cp 2 at time 50 minutes (3000 seconds) then stop the
    # maintenance infusions
    p1.add_target(start=3000,
                  target=2,
                  duration=10,
                  effect=False,
                  maintenance_infusions=False)

    # Simulated plasma and effect site concentrations
    concentrations = p1.run()

    print("\nConcentrations"
          "\n=============="
          "\nColumn 0: Time (seconds)"
          "\nColumn 1: Plasma concentration"
          "\nColumn 2: Effect site concentration")

    print(concentrations)

    # Infusion List
    print("\nInfusion list"
          "\n============="
          "\nColumn 0: Start time (seconds)"
          "\nColumn 1: Dose (mg) per second"
          "\nColumn 2: Duration (seconds)"
          "\nColumn 3: End time (seconds)\n")

    print(p1.infusion_list)

    # Rates list
    rates = p1.generate_rates_array()

    print("\nRates"
          "\n====="
          "\nColumn 0: Time (seconds)"
          "\nColumn 1: Rate (ml/hr)\n")

    print(rates)

   
