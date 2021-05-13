
"""
opentiva.pump
=============

This module contains the class for simulating the pump.
"""

import warnings
import numpy as np
import scipy.optimize as optimize

import opentiva.pkpd as pkpd


class Pump:
    """Pump class simulates a target controlled infusion pump for a drug's
    pharmacokinetic/ pharmacodynamic model.

    Parameters
    ----------
    model
        opentiva pharmacokinetic/ pharmacodynamic model.
        e.g. opentiva.propofol.MarshDiprifusor(...)
    drug_concentration
        concentration of infusion drug
        (using units outlined in model; model.concentration_unit)
    end_time
        duration of simulation in seconds
    cp_limit
        multiplied to the an effect site target to get maximum increase
        in plasma concentration in 'revised' effect site targetting
        value provided here is default for the targets.
    cp_limit_duration
        time to achieve cp_limit in seconds
        value provided here is default for the targets
    maintenance_infusion_duration
        time in seconds of the duration of each maintenance
        infusion
    maintenance_infusion_multiplier
        duration of each subsequent maintenance infusion is multiplied
        by this, extending the time between each maintenance infusion
        calculation
    max_infusion_rate
        ml/hr limit on infusion rate of pump
    bolus_time
        time in seconds below which infusions are considered as a 'bolus'
        i.e. not effected by the max_infusion_rate

    Attributes
    ----------
    target_concentrations
        2d array of targets with each row containing:

           #. start:  start time of target in seconds
           #. target:  target of plasma or effect concentration
           #. duration:  time to achieve target level in seconds
           #. effect:  true for effect site targetting or false for plasma
              site targetting
           #. cp_limit:  multiplied to the target to get maximum increase in
              plasma concentration in 'revised' effect site targetting
           #. cp_limit_duration:  time to achieve cp_limit in seconds
           #. ce_bolus_only:  true will target the effect using the
              'original' method / bolus only method
           #. maintenance_infusions: true will calculate the infusions between
              the time the target is reached and the next target to maintain
              steady state

    infusion_list
        2d array of infusions calculated from provided targets
        with each row containing:
           #. start time of infusion in seconds
           #. dose of infusion over 1 second
           #. duration of infusion in seconds
           #. end time of infusion in seconds
    user_infusion_list
        2d array of infusions entered by user
        with each row containing:
           #. start time of infusion in seconds
           #. dose of infusion over 1 second
           #. duration of infusion in seconds
           #. end time of infusion in seconds
    """

    def __init__(self, model, drug_concentration: int, end_time: int,
                 cp_limit: float = 1.2,
                 cp_limit_duration: int = 10,
                 maintenance_infusion_duration: int = 300,
                 maintenance_infusion_multiplier: int = 2,
                 max_infusion_rate: int = 1200,
                 bolus_time: int = 20):

        # Instance variables
        self.model = model
        self.drug_concentration = drug_concentration
        self.end_time = end_time
        self.cp_limit = cp_limit
        self.cp_limit_duration = cp_limit_duration
        self.maintenance_infusion_multiplier = maintenance_infusion_multiplier
        self.maintenance_infusion_duration = maintenance_infusion_duration
        self.max_infusion_rate = max_infusion_rate
        self.bolus_time = bolus_time

        # Arrays used in class
        self.target_concentrations = np.empty((0, 9))
        self.infusion_list = np.empty((0, 4))
        self.user_infusion_list = np.empty((0, 4))
        self.mi_delta = np.zeros(1)

        # Pharmacokinetic / pharmacokdynamic model
        self.pkpd_model = pkpd.PkPdModel(self.model)

    # Validation instance variables

    @property
    def drug_concentration(self):
        return self._drug_concentration

    @drug_concentration.setter
    def drug_concentration(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError("Expected: float or int")
        if value <= 0:
            raise ValueError("Must be greater than 0")
        self._drug_concentration = value

    @property
    def end_time(self):
        return self._end_time

    @end_time.setter
    def end_time(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError("Expected: float or int")
        if value <= 0:
            raise ValueError("Must be greater than 0")
        self._end_time = value

    @property
    def cp_limit(self):
        return self._cp_limit

    @cp_limit.setter
    def cp_limit(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError("Expected: float or int")
        if value <= 0:
            raise ValueError("Must be greater than 0")
        self._cp_limit = value

    @property
    def cp_limit_duration(self):
        return self._cp_limit_duration

    @cp_limit_duration.setter
    def cp_limit_duration(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError("Expected: float or int")
        if value <= 0:
            raise ValueError("Must be greater than 0")
        self._cp_limit_duration = value

    @property
    def maintenance_infusion_multiplier(self):
        return self._maintenance_infusion_multiplier

    @maintenance_infusion_multiplier.setter
    def maintenance_infusion_multiplier(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError("Expected: float or int")
        if value <= 0:
            raise ValueError("Must be greater than 0")
        self._maintenance_infusion_multiplier = value

    @property
    def maintenance_infusion_duration(self):
        return self._maintenance_infusion_duration

    @maintenance_infusion_duration.setter
    def maintenance_infusion_duration(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError("Expected: float or int")
        if value <= 0:
            raise ValueError("Must be greater than 0")
        self._maintenance_infusion_duration = value

    @property
    def max_infusion_rate(self):
        return self._max_infusion_rate

    @max_infusion_rate.setter
    def max_infusion_rate(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError("Expected: float or int")
        if value <= 0:
            raise ValueError("Must be greater than 0")
        self._max_infusion_rate = value

    @property
    def bolus_time(self):
        return self._bolus_time

    @bolus_time.setter
    def bolus_time(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError("Expected: float or int")
        if value <= 0:
            raise ValueError("Must be greater than 0")
        self._bolus_time = value

    # Class methods

    def _add_infusion(self, start: int, dose: float, duration: int) -> None:
        """Adds an infusion to infusion_list
        """
        end = start + duration
        inf = np.array([[start, dose, duration, end]])
        self.infusion_list = np.vstack((self.infusion_list, inf))

    def add_infusion(self, start: int, dose: float, duration: int) -> None:
        """ Adds a user defined infusion to the infusion list

        Parameters
        ----------
        start
            start time of infusion in seconds
        dose
            dose per second (e.g. mg / second) of infusion
        duration
            total time of infusion in seconds

        Returns
        -------
        None
            Adds infusion to user infusion list
        """
        end = start + duration
        inf = np.array([[start, dose, duration, end]])
        self.user_infusion_list = np.vstack((self.user_infusion_list, inf))

    def add_target(self, start: int, target: float, duration: int,
                   effect: bool, cp_limit: float = 0,
                   cp_limit_duration: float = 0, ce_bolus_only: bool = True,
                   maintenance_infusions: bool = True) -> None:
        """Takes parameters for a target and adds it to the
        target_concentrations array

        Parameters
        ----------
        start
            start time of target in seconds
        target
            target of plasma or effect concentration
        duration
            time to achieve target level in seconds
        effect
            true for effect site targetting or false for plasma site targetting
        cp_limit
            multiplied to the target to get maximum increase in plasma
            concentration in 'revised' effect site targetting
        cp_limit_duration
            time to achieve cp_limit in seconds
        ce_bolus_only
            true will target the effect using the 'original' method / bolus
            only method
        maintenance_infusions
            true will calculate the infusions between the time the target is
            reached and the next target to maintain steady state

        Returns
        -------
        None
            Adds above target parameters to target_concentrations array
        """

        if cp_limit == 0:
            cp_limit = self.cp_limit
        if cp_limit_duration == 0:
            cp_limit_duration = self.cp_limit_duration

        tc_arr = self.target_concentrations
        tc_len = self.target_concentrations.shape[0]

        target_v = np.array([start, target, duration, 0, effect, cp_limit,
                            cp_limit_duration, ce_bolus_only,
                            maintenance_infusions])
        tc_arr = np.vstack((tc_arr, target_v))
        tc_arr_sorted = tc_arr[tc_arr[:, 0].argsort()]

        for n in range(tc_len):
            # Iterate through tc_arr_sorted and change end times to the next
            # infusions start time minus one
            tc_arr_sorted[n, 3] = tc_arr_sorted[n + 1, 0] - 1

        # Final target's end time to match end time of simulation
        tc_arr_sorted[-1, 3] = self.end_time
        self.target_concentrations = tc_arr_sorted

    def generate_infusions(self) -> None:
        """ Generates the infusions required to achieve the targets in the
        target concentrations array and stores in infusion_list
        """

        self.infusion_list = np.empty((0, 4))
        self.mi_delta = np.zeros(1)
        tc_arr = self.target_concentrations
        tc_len = self.target_concentrations.shape[0]

        for n in range(tc_len):
            start = tc_arr[n, 0]
            target = tc_arr[n, 1]
            duration = tc_arr[n, 2]
            end = tc_arr[n, 3]
            effect = tc_arr[n, 4]
            cp_limit = tc_arr[n, 5]
            cp_limit_duration = tc_arr[n, 6]
            ce_bolus_only = tc_arr[n, 7]
            maintenance_infusions = tc_arr[n, 8]

            if n == 0:
                # Initial target
                self._concentration_increase(start, target, duration, end,
                                             effect, cp_limit,
                                             cp_limit_duration,
                                             ce_bolus_only,
                                             maintenance_infusions)
            else:
                # Other targets determine if increase/ decrease in
                # concentration
                c_delta = np.diff((tc_arr[:, 1]))[n - 1]
                if c_delta > 0:
                    self._concentration_increase(start, target, duration, end,
                                                 effect, cp_limit,
                                                 cp_limit_duration,
                                                 ce_bolus_only,
                                                 maintenance_infusions)
                else:
                    # concentration decrease
                    self._concentration_decrease(start, target, duration, end,
                                                 effect, maintenance_infusions)

        # Add user specified infusions
        self.infusion_list = np.vstack((self.infusion_list,
                                        self.user_infusion_list))

    def _concentration_increase(self, start: int, target: float,
                                duration: int, end_target: int, effect: bool,
                                cp_limit: float, cp_limit_duration: int,
                                ce_bolus_only: bool,
                                maintenance_infusions: bool) -> None:
        """Method handles generating the infusions relating to a target
        concentration increase
        """

        end = start + duration

        # Get previous concentration
        if start == 0:
            previous_cp = 0
        else:
            previous_cp, _ = self.get_concentration(start - 1)

        if effect:
            break_count = 0
            while True:
                if ce_bolus_only:
                    root = optimize.newton(self.pkpd_model.ce_cplimit_minimise,
                                           x0=1, x1=10,
                                           args=(self.infusion_list,
                                                 target,
                                                 cp_limit_duration,
                                                 start,
                                                 duration,
                                                 self.drug_concentration,
                                                 self.max_infusion_rate,
                                                 self.bolus_time))

                    cp_limit = root

                    # Update tc array with calculated Cp Limit
                    for tc_v in self.target_concentrations:
                        if tc_v[0] == start:
                            tc_v[5] = cp_limit
                            tc_v[6] = cp_limit_duration

                inf_out, target_time = self.pkpd_model.ce_dose(self.infusion_list,
                                                               target,
                                                               cp_limit,
                                                               cp_limit_duration,
                                                               start,
                                                               duration,
                                                               self.drug_concentration,
                                                               self.max_infusion_rate,
                                                               self.bolus_time)

                # If specified duration is greater than time to target
                # default to ce_bolus_only and prolong cp_limit_duration
                if target_time < end:
                    cp_limit_duration = duration
                    ce_bolus_only = True
                elif break_count > 2:
                    break
                else:
                    break

                break_count += 1

            self.infusion_list = inf_out
            end = target_time
        else:
            c = target - previous_cp
            rate = self.max_infusion_rate + 1

            while True:
                dose = c / self.pkpd_model.integral_exp_decline(0, duration)
                rate = (dose / self.drug_concentration) * 60 * 60

                if duration <= self.bolus_time:
                    break
                elif self.max_infusion_rate == -1:
                    break
                elif rate <= self.max_infusion_rate:
                    break
                else:
                    duration += 1

                if duration > 100:
                    break

            self._add_infusion(start, dose, duration)
            end = start + duration

        if maintenance_infusions:
            # Generate maintenance infusion after loading dose complete
            tc_v = np.array([end, target, end_target])
            self._generate_maintenance_infusions(tc_v)
        else:
            # Add zero dose infusion if no maintenance infusions
            self._add_infusion(end, 0, end_target)

    def _concentration_decrease(self, start: int, target: float, duration: int,
                                end_target: int, effect: bool,
                                maintenance_infusions: bool) -> None:
        """Method handles generating the infusions relating to a target
        concentration decrease
        """

        if effect:
            dec_time = self.pkpd_model.effect_decrement_time(start, target,
                                                             self.infusion_list)
        else:
            dec_time = self.pkpd_model.plasma_decrement_time(start, target,
                                                             self.infusion_list)

        if duration < dec_time:
            # If duration is less than time for natural expontential decline
            # then end time of 0 infusions matches decrement time
            dose = 0
            self._add_infusion(start, dose, dec_time)
            end = start + dec_time
        else:
            # If duration is greater than natural expontential decline add
            # infusion to offset to decline to reach specified duration

            dec_cp, dec_ce = self.get_concentration(start + duration)

            if effect:
                delta_c = target - dec_ce
            else:
                delta_c = target - dec_cp

            dose = delta_c / self.pkpd_model.integral_exp_decline(0, duration)
            self._add_infusion(start, dose, duration)
            end = start + duration

        if maintenance_infusions:
            # Generate maintenance infusion after loading dose complete
            tc_v = np.array([end, target, end_target])
            self._generate_maintenance_infusions(tc_v)
        else:
            # Add zero dose infusion if no maintenance infusions
            self._add_infusion(end, 0, end_target)

    def _generate_maintenance_infusions(self, target_v) -> None:
        """Method handles generating the infusions required to offset the
        clearence and elimation loses to maintain a steady state
        """

        inf_out = self.pkpd_model.maintenance_infusion_list(target_v,
                                                            self.infusion_list,
                                                            self.maintenance_infusion_duration,
                                                            self.maintenance_infusion_multiplier,
                                                            self.drug_concentration,
                                                            self.max_infusion_rate)
        self.infusion_list = inf_out

    def decrement_time(self, start: int, effect: bool, target: float) -> int:
        """ Returns time in seconds to reach a target from a start time
        (seconds) if all infusions were stopped

        Parameters
        ----------
        start
            start time to calculate decrement time from in seconds
        effect
            true for effect site decrement time or false for plasma site
            decrement time
        target
            target concentration of plasma or effect concentration to
            decrement to

        Returns
        -------
        int
            time in seconds to reach target concentration
        """

        decrement_time = 0

        if effect:
            decrement_time = self.pkpd_model.effect_decrement_time(start,
                                                                  target,
                                                                  self.infusion_list)
        elif not effect:
            decrement_time = self.pkpd_model.plasma_decrement_time(start,
                                                                  target,
                                                                  self.infusion_list)
        return decrement_time

    def generate_rates_array(self) -> np.ndarray:
        """ Method turns the infusion_list array into an output ml/hr array

        This does not include user defined infusions

        Parameters
        ----------
        None

        Returns
        -------
        np.ndarray
            column 0: time of rate change (seconds)
            column 1: rate (ml/hr)
        """

        if (self.infusion_list.size - self.user_infusion_list.size) == 0:
            # if infusion_list contains only user defined infusions return None
            return np.array([])

        x_max = self.infusion_list.shape[0] - self.user_infusion_list.shape[0]
        #  Ignore used defined infusions
        rates_arr = np.zeros((x_max + 1, 2))

        for x in range(x_max):
            dose = self.infusion_list[x, 1]
            start = self.infusion_list[x, 0]
            rate = (dose / self.drug_concentration) * 60 * 60
            # Dose is 'dose per second' in infusion_list

            rates_arr[x, 0] = start
            rates_arr[x, 1] = rate

        rates_arr[-1, 0] = self.end_time
        rates_arr[-1, 1] = rates_arr[-2, 1]

        return rates_arr

    def generate_dose_weight_array(self, interval: str = 'min') -> np.ndarray:
        """ Method turns the infusion_list array into dose/weight (if infusion
        time below bolus time) or dose/weight/time array if not

        Parameters
        ----------
        interval
            'min' or 'hr', gives time interval for
            dose/weight/time output. Default 'min'

        Returns
        -------
        np.ndarray
            column 0: time of rate change (seconds);
            column 1: boluses dose expressed as dose/weight (e.g. mg/kg) or
            dose/weight/time (e.g. mg/kg/hr);
            column 2: true if bolus dose; dose/weight (e.g. mg/kg)
            false if infusion; dose/weight/time (e.g. mg/kg/hr)
        """

        x_max = self.infusion_list.shape[0]
        dose_weight_arr = np.zeros((x_max + 1, 3))
        weight = self.model.weight

        if interval == 'min':
            time_interval = 60
        elif interval == 'hr':
            time_interval = 60 * 60
        else:
            time_interval = 60
            warnings.warn("Warning: func generate_dose_weight_array. Use 'min'\
                           or 'hr' for interval. Defaulting to 'min'")

        for x in range(x_max):
            start = self.infusion_list[x, 0]
            dose = self.infusion_list[x, 1]
            duration = self.infusion_list[x, 2]

            if duration <= self.bolus_time and \
                    duration < self.maintenance_infusion_duration:
                dose_weight = (dose * duration) / weight
                dose_bolus = True
            else:
                dose_weight = dose / weight * time_interval
                dose_bolus = False

            dose_weight_arr[x, 0] = start
            dose_weight_arr[x, 1] = dose_weight
            dose_weight_arr[x, 2] = dose_bolus

        dose_weight_arr[-1, 0] = self.end_time
        dose_weight_arr[-1, 1] = dose_weight_arr[-2, 1]
        dose_weight_arr[-1, 2] = dose_weight_arr[-2, 2]

        return dose_weight_arr

    def generate_targets_array(self) -> np.ndarray:
        """ Method turns the target_concentrations array to a targets array

        Parameters
        ----------
        None

        Returns
        -------
        np.ndarray
            column 0: time of target change (seconds);
            column 1: target concentration
        """

        x_max = self.target_concentrations.shape[0]
        targets_arr = np.zeros((x_max + 1, 2))

        for x in range(x_max):
            start = self.target_concentrations[x, 0]
            target = self.target_concentrations[x, 1]

            targets_arr[x, 0] = start
            targets_arr[x, 1] = target

        targets_arr[-1, 0] = self.end_time
        targets_arr[-1, 1] = targets_arr[-2, 1]

        return targets_arr

    def get_concentration(self, time: int):
        """ Method returns the plasma and effect site concentration at point
        in time (seconds)

        Parameters
        ----------
        time
            time to get concentrations at in seconds

        Returns
        -------
        tuple[float, float]
            plasma concentration and effect site concentration
        """
        cp = self.pkpd_model.cp_over_time(self.infusion_list, 0, time)
        ce = self.pkpd_model.ce_over_time(cp)
        return cp[-1, 1], ce[-1][0]

    def run(self) -> np.ndarray:
        """ Method returns plasma and effect concentrations over the simulation
        time period

        Parameters
        ----------
        None

        Returns
        -------
        np.ndarray
            column 0: time of concentration in seconds
            column 1: plasma site concentration
            column 2: effect site concentration
        """

        self.generate_infusions()

        cp_arr = self.pkpd_model.cp_over_time(self.infusion_list, 0,
                                              int(self.end_time))
        ce_arr = self.pkpd_model.ce_over_time(cp_arr)

        concentrations = np.hstack((cp_arr, ce_arr))
        return concentrations
