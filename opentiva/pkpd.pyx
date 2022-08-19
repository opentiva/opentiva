import warnings
import numpy as np
import scipy.integrate as integrate
import scipy.optimize as optimize

from libc.math cimport sqrt
from libc.math cimport acos
from libc.math cimport cos
from libc.math cimport pi
from libc.math cimport exp
from libc.math cimport log


cdef class PkPdModel:
    """PkPdModel class is used to model the pharmacokinetic and pharmacodyamics
    mathematical models

    Parameters
    ----------
    model
        opentiva pharmacokinetic/ pharmacodynamic model.
        e.g. opentiva.propofol.MarshDiprifusor(...)

    """

    cdef double alpha, beta, gamma, A, B, C, v1
    cdef double k10, k12, k13, k21, k31, k20, ke0

    def __init__(self, model):

        # Convert model's rate constants to seconds
        self.v1 = model.v1
        self.k10 = model.k10 / 60
        self.ke0 = model.ke0 / 60
        
        cdef int compartments = model.compartments

        if compartments == 3:
            self.k13 = model.k13 / 60
            self.k31 = model.k31 / 60
            self.k12 = model.k12 / 60
            self.k21 = model.k21 / 60
            self.three_compartment()
        elif compartments == 2:
            self.k12 = model.k12 / 60
            self.k21 = model.k21 / 60

            if hasattr(model, 'k20'):
                self.k20 = model.k20 / 60
            else:
                self.k20 = 0

            self.two_compartment()
        elif compartments == 1:
            self.one_compartment()
        else:
           raise ValueError("Compartment variables must be 1, 2 or 3")

    cdef three_compartment(self):
        """Adds the three compartment model variables to the class'
        instance variables

        Reference: http://www.pfim.biostat.fr/PFIM_PKPD_library.pdf
        """
        cdef double a0, a1, a2, p, q, r1, r2, theta

        # Three compartment model with linear elimination variables

        a0 = self.k10 * self.k21 * self.k31
        a1 = (self.k10 * self.k31) + (self.k21 * self.k31) \
            + (self.k21 * self.k13) + (self.k10 * self.k21) + \
            (self.k31 * self.k12)
        a2 = self.k10 + self.k12 + self.k13 + self.k21 + self.k31

        p = a1 - (a2 ** 2 / 3)
        q = ((2 * a2 ** 3) / 27) - (a1 * a2 / 3) + a0

        r1 = sqrt(-(p ** 3 / 27))
        r2 = 2 * r1 ** (1 / 3.0)

        theta = acos(-(q / (2 * r1))) / 3

        self.alpha = -(cos(theta) * r2 - a2 / 3)
        self.beta = -(cos(theta + (2 * pi) / 3) *
                      r2 - a2 / 3)
        self.gamma = -(cos(theta + (4 * pi) / 3) *
                       r2 - a2 / 3)

        self.A = (1 / self.v1) * \
            ((self.k21 - self.alpha) / (self.alpha - self.beta)) * \
            ((self.k31 - self.alpha) / (self.alpha - self.gamma))
        self.B = (1 / self.v1) * \
            ((self.k21 - self.beta) / (self.beta - self.alpha)) * \
            ((self.k31 - self.beta) / (self.beta - self.gamma))
        self.C = (1 / self.v1) * \
            ((self.k21 - self.gamma) / (self.gamma - self.beta)) * \
            ((self.k31 - self.gamma) / (self.gamma - self.alpha))

    cdef two_compartment(self):
        """Adds the two compartment model variables to the class'
        instance variables with optional k20 elimination

        Reference: http://www.pfim.biostat.fr/PFIM_PKPD_library.pdf
        """
        # Two compartment model with linear elimination variables and 
        # optional k20 elimination
        a1 = (self.k21 * self.k10) + (self.k12 * self.k20) + \
             (self.k10 * self.k20)
        a2 = self.k12 + self.k21 + self.k10 + self.k20

        self.beta = 0.5 * (a2 - sqrt(a2 ** 2 - (4 * a1)))
        self.alpha = a1 / self.beta
        self.gamma = 1  # arbitrary set gamma as it will be ignored as C = 0

        self.A = 1 / self.v1 * \
                ((self.alpha - self.k21 - self.k20) / (self.alpha - self.beta))
        self.B = 1 / self.v1 * \
                ((self.beta - self.k21 - self.k20) / (self.beta - self.alpha))
        self.C = 0

    cdef one_compartment(self):
        """Adds the one compartment model variables to the class'
        instance variables

        Reference: http://www.pfim.biostat.fr/PFIM_PKPD_library.pdf
        """
        # One compartment model with linear elimination variables
        self.alpha = self.k10

        self.beta = 1  # arbitrary set beta as it will be ignored as B = 0
        self.gamma = 1  # arbitrary set gamma as it will be ignored as C = 0

        self.A = 1 / self.v1
        self.B = 0
        self.C = 0

    cpdef double integrand_exp_decline(self, double time):
        """ Method returns value of the three compartment exponential decline
        function at a point in time
        """
        cdef double f = (self.A * exp(-self.alpha * time) + \
                         self.B * exp(-self.beta * time) + \
                         self.C * exp(-self.gamma * time))
        return f


    cpdef double integral_exp_decline(self, double x_min, double x_max):
        """ Method integrates the exponential decline function over time
        """
        cdef tuple i = integrate.quad(self.integrand_exp_decline, x_min, x_max)
        return i[0]


    cdef double cp_increment(self, double dose, int elapsed):
        """Returns the increase in plasma concentration from an
        infusion runnning over a time period

        Parameters
        ----------
        dose
            dose of drug over 1 second
        elapsed
            time in seconds since infusion started

        Returns
        -------
        double
            plasma concentration of increment
        """

        cdef double cp_inc  = dose * (
            (self.A / self.alpha) * (1 - exp(-self.alpha * elapsed)) +
            (self.B / self.beta) * (1 - exp(-self.beta * elapsed)) +
            (self.C / self.gamma) * (1 - exp(-self.gamma * elapsed)))

        return cp_inc


    cdef double cp_decrement(self, double dose, int duration, int elapsed):
        """ Returns the plasma concentration after an infusion has been
        stopped at a point in time representing the clearance and elimation
        loses

        Parameters
        ----------
        dose
            dose of drug over 1 second
        duration
            time in seconds of infusion duration
        elapsed
            time in seconds since infusion ended

        Returns
        -------
        double
            plasma concentration of decrement
         """

        cdef double cp_dec, a, b, c

        a = ((self.A / self.alpha) * (1 - exp(-self.alpha * duration)) * \
                (exp(-self.alpha * elapsed)))

        b = ((self.B / self.beta) * (1 - exp(-self.beta * duration)) * \
                (exp(-self.beta * elapsed)))

        c = ((self.C / self.gamma) * (1 - exp(-self.gamma * duration)) * \
                (exp(-self.gamma * elapsed)))

        cp_dec = dose * (a + b + c)

        return cp_dec


    cpdef double calculate_cp(self, double [:, :] infusion_list, int time):
        """Takes an array of infusions and returns the plasma concentration
        at a point in time

        Parameters
        ----------
        infusion_list
            2d array of infusions with each row containing:
            [start time of infusion in seconds,
            dose of infusion over 1 second,
            duration of infusion in seconds,
            end time of infusion in seconds]
        time
            time in seconds to calculate plasma concentration at

        Returns
        -------
        double
            plasma concentration at a point in time
        """

        cdef double dose, cp, target
        cdef int start, duration, end, elapsed, diff

        cdef Py_ssize_t x_max = int(infusion_list.shape[0])
        cdef Py_ssize_t x

        cp = 0

        for x in range(x_max):
            start = int(infusion_list[x, 0])
            dose = infusion_list[x, 1]
            duration = int(infusion_list[x, 2])
            end = int(infusion_list[x, 3])
            elapsed = time - start  # Time since infusion started running
            diff = time - end  # Time since infusion stopped

            if (time <= end) and (time >= start):
                cp += self.cp_increment(dose, elapsed)
            elif time > end:
                cp += self.cp_decrement(dose, duration, diff)

        return cp


    cpdef cp_over_time(self, double[:, :] infusion_list, int start, int end):
        """Takes an array of infusions and returns the plasma concentrations
        over a time range

        Parameters
        ----------
        infusion_list
            2d array of infusions with each row containing:
            [start time of infusion in seconds,
            dose of infusion over 1 second,
            duration of infusion in seconds,
            end time of infusion in seconds]
        start
            time in seconds to calculate plasma concentration from
        end
            time in seconds to calculate plasma concentration til

        Returns
        -------
        np.ndarray
            2d array of plasma concentrations over time with each row
            containing:
            [time of plasma concentration in seconds,
            plasma concentration]
        """

        cdef int x = 0, delta
        cdef double cp
        cdef Py_ssize_t t

        delta = end - start
        cp_arr = np.empty((delta, 2), dtype=np.float64)

        for t in range(start, end):
            cp = self.calculate_cp(infusion_list, t)

            cp_arr[x, 0] = t
            cp_arr[x, 1] = cp
            x += 1

        return cp_arr


    cdef double calculate_ce(self, double current_cp, double previous_cp,
                             double previous_ce):
        """Returns the effect site concentration at a point in time

        Parameters
        ----------
        current_cp
            plasma concentration at time
        previous_cp
            plasma concentration at time - 1
        previous_ce
            effect site concentration at time - 1

        Returns
        -------
        double
            effect site concentration at a point in time
        """

        cdef double current_ce, delta_cp
        cdef double delta = 0

        delta_cp = current_cp - previous_cp

        if previous_cp == 0:
            return 0  # avoid divide by zero error

        if delta_cp > 0:
            slope = delta_cp
            delta =  (1 * slope + (self.ke0 * previous_cp - slope)) * \
                     (1 - exp(-self.ke0 * 1)) / self.ke0

        elif delta_cp <= 0:
            slope = log(current_cp) - log(previous_cp)
            delta =  previous_cp * self.ke0 / (self.ke0 + slope) * \
                     (exp(1 * slope) - exp(-self.ke0 * 1))

        current_ce = previous_ce * exp(-self.ke0) + delta

        return current_ce


    cpdef ce_over_time(self, double[:, :] cp_arr):
        """Takes an array of plasma concentrations starting from time 0
        and returns the effect site concentrations over that time range

        Parameters
        ----------
        cp_arr
            2d array of plasma concentrations over time with each row
            containing (from cp_over_time function starting at time 0):
            [time of plasma concentration in seconds,
            plasma concentration]

        Returns
        -------
        np.ndarray
            1d array of effect site concentration over time
        """

        cdef double current_cp, previous_cp, delta_cp, current_ce, previous_ce
        cdef Py_ssize_t x_max = int(cp_arr.shape[0])
        cdef Py_ssize_t x

        ce = np.zeros((x_max, 1), dtype=np.float64)

        previous_ce = 0
        current_ce = 0

        for x in range(x_max):

            if x == 0:
                continue  # skip first cp

            previous_cp = cp_arr[x - 1, 1]
            current_cp = cp_arr[x, 1]

            current_ce = self.calculate_ce(current_cp, previous_cp,
                                           previous_ce)

            previous_ce = current_ce

            ce[x] = current_ce

        return ce


    cpdef ce_dose(self, double [:, :] infusion_list, double target,
                  double limit, int duration_b, int start_b, int duration_ce,
                  double drug_concentration, int max_infusion_rate,
                  int bolus_time):
        """Returns the infusions require to reach a target effect site
        concentration

        Parameters
        ----------
        infusion_list
            2d array of infusions with each row containing:
            [start time of infusion in seconds,
            dose of infusion over 1 second,
            duration of infusion in seconds,
            end time of infusion in seconds]
        target
            effect site concentration to reach
        limit
            multiplied by the target to give the max plasma concentration
            during the targetting
        duration_b
            time in seconds to reach initial limit plasma target concentration
        start_b
            time in seconds of the start of the effect site targetting
        duration_ce
            time in seconds to reach effect site target
        drug_concentration
            concentration of infusion drug
        max_infusion_rate
            ml/hr limit on infusion rate of pump
        bolus_time
            time in seconds below which infusions are considered as a 'bolus'
            i.e. not effected by the max_infusion_rate

        Returns
        -------
        np.ndarray
            2d array of infusions with each row containing:
            [start time of infusion in seconds,
            dose of infusion over 1 second,
            duration of infusion in seconds,
            end time of infusion in seconds]
        int
            time in seconds that effect site concentration reached
        """

        cdef double target_limit, root, rate
        cdef double previous_cp, delta_cp, dose_cp, dose_mi
        cdef int start_mi, duration_mi, end_mi, end_b
        cdef int target_time
        inf_out = infusion_list

        # Get dose to increment to max Cp limit over limit duration
        if start_b == 0:
            previous_cp = 0
        else:
            previous_cp = self.calculate_cp(inf_out, start_b)

        target_limit = target * limit
        delta_cp = target_limit - previous_cp

        while True:  # Extend bolus dose to max infusion rate
            dose_cp = delta_cp / self.integral_exp_decline(0, duration_b)
            rate = (dose_cp / drug_concentration) * 60 * 60

            if duration_b <= bolus_time:
                break
            elif max_infusion_rate == -1:
                break
            elif rate <= max_infusion_rate:
                break
            else:
                duration_b += 1

        # Add starting bolus dose to array
        end_b = start_b + duration_b
        inf_v = np.array([start_b, dose_cp, duration_b, end_b], dtype=np.float64)
        inf_out = np.vstack((inf_out, inf_v))

        # Use newton secant to find duration required for Ce to reach target
        start_mi = end_b

        try:
            root = optimize.newton(self.ce_duration_minimise,
                                   x0 = 1, x1 = duration_b * 2, tol=1,
                                   args=(inf_out, target, limit, start_mi))
        except (RuntimeError, OverflowError):
            warnings.warn("Failed to converge on infusion time.", RuntimeWarning)
            end_mi = end_b
        else:
            duration_mi = int(root)

            # Stop negative durations
            if duration_mi < 0:
                duration_mi = 0

            dose_mi = self.maintenance_infusion(inf_out, target_limit,
                                                start_mi, duration_mi)
            end_mi = start_mi + duration_mi

            # If rate above max rate match match max infusion rate
            rate = (dose_mi / drug_concentration) * 60 * 60
            if rate > max_infusion_rate:
                dose_mi = rate / (60 * 60) * drug_concentration

            inf_v = np.array([start_mi, dose_mi, duration_mi, end_mi],
                              dtype=np.float64)

            # Add infusion if duration_mi > 0
            if duration_mi:
                inf_out = np.vstack((inf_out, inf_v))

        # Find time at which Ce reaches target
        target_time = end_mi
        cp = target + limit
        while cp >= target:
            cp = self.calculate_cp(inf_out, target_time)
            target_time += 1

        # Add zero infusion til ce reached if maintenance infusion required
        inf_0 = np.array([end_mi, 0, (target_time - end_mi), target_time])
        inf_out = np.vstack((inf_out, inf_0))

        return inf_out, target_time


    cpdef double ce_duration_minimise(self, int duration, 
                                      double [:, :] infusion_list, 
                                      double target, double limit, int start):
        """ Minimisation function to calculate duration of Tinf

        Method when minimised by changing the duration variable
        will give the duration of the Tinf phase used to to reach an
        effect site target using the method described by
        Van Poucke et al (2004, PMID: 15536889 DOI: 10.1109/TBME.2004.827935)

        Parameters
        ----------
        duration
            duration in seconds of Tcoast phase of the revised effect site
            targeting
        infusion_list
            2d array of infusions with each row containing:
            [start time of infusion in seconds,
            dose of infusion over 1 second,
            duration of infusion in seconds,
            end time of infusion in seconds]
        target
            effect site concentration to reach
        limit
            multiplied by target to give the max plasma concentration
            during the targetting
        start
            time in seconds of the start of the effect site targetting

        Returns
        -------
        double
            minimisation target

        """

        cdef double target_limit, current_ce, previous_ce, delta_ce
        cdef double current_cp, previous_cp
        cdef int t = 1, end

        target_limit = target * limit
        previous_ce = 0
        current_cp = 0
        end = start + duration

        dose = self.maintenance_infusion(infusion_list, target_limit,
                                         start, duration)

        inf_v = np.array([start, dose, duration, end], dtype=np.float64)
        inf_tmp = np.vstack((infusion_list, inf_v))

        while True:
            previous_cp = current_cp
            current_cp = self.calculate_cp(inf_tmp, t)

            current_ce = self.calculate_ce(current_cp, previous_cp,
                                           previous_ce)

            delta_ce = previous_ce - current_ce
            previous_ce = current_ce

            if delta_ce >= 0 and  t > end:
                break

            t += 1

        return target - current_ce


    cpdef ce_cplimit_minimise(self, double limit, double [:, :] infusion_list,
                  double target, int duration_b, int start_b,
                  int duration_ce, double drug_concentration,
                  int max_infusion_rate, int bolus_time):
        """Minimisation function to calculate limit value for original targeting

        Method when minimised by changing the limit variable will give
        the limit value that when multiplied to the effect target concentration
        gives the required maxmimum plasma concentration t to achieve an
        effect target using the original Shafer/Gregg method
        (1992 PMID: 1629794 DOI: 10.1007/BF01070999)

        Parameters
        ----------
        limit
            multiplied by the target to give the max plasma concentration
            during the targetting
        infusion_list
            2d array of infusions with each row containing:
            [start time of infusion in seconds,
            dose of infusion over 1 second,
            duration of infusion in seconds,
            end time of infusion in seconds]
        target
            effect site concentration to reach
        duration_b
            time in seconds to reach initial limit plasma target concentration
        start_b
            time in seconds of the start of the effect site targetting
        duration_ce
            time in seconds to reach effect site target
        drug_concentration
            concentration of infusion drug
        max_infusion_rate
            ml/hr limit on infusion rate of pump
        bolus_time
            time in seconds below which infusions are considered as a 'bolus'
            i.e. not effected by the max_infusion_rate

        Returns
        -------
        double
            minimisation target
        """

        cdef double current_cp, previous_cp, delta_cp
        cdef double current_ce, previous_ce, delta_ce
        cdef int start_mi, duration_mi, end_mi, end_b
        cdef int t = 1

        inf_tmp = infusion_list

        # Get dose to increment to max Cp limit over limit duration
        if start_b == 0:
            previous_cp = 0
        else:
            previous_cp = self.calculate_cp(inf_tmp, start_b)

        target_limit = target * limit
        delta_cp = target_limit - previous_cp

        while True:  # Extend bolus dose to max infusion rate
            dose_cp = delta_cp / self.integral_exp_decline(0, duration_b)
            rate = (dose_cp / drug_concentration) * 60 * 60

            if duration_b <= bolus_time:
                break
            elif max_infusion_rate == -1:
                break
            elif rate <= max_infusion_rate:
                break
            else:
                duration_b += 1

        # Add bolus dose to array
        end_b = start_b + duration_b
        inf_v = np.array([start_b, dose_cp, duration_b, end_b], 
                          dtype=np.float64)
        inf_tmp = np.vstack((inf_tmp, inf_v))

        previous_ce = 0
        current_cp = 0
        while True:
            previous_cp = current_cp
            current_cp = self.calculate_cp(inf_tmp, t)

            current_ce = self.calculate_ce(current_cp, previous_cp,
                                           previous_ce)

            delta_ce = previous_ce - current_ce
            previous_ce = current_ce

            if delta_ce > 0 and  t >= end_b:
                break

            t += 1

        return target - current_ce


    cpdef maintenance_infusion(self, double [:, :] infusion_list, double target,
                               int time, int duration):
        """Returns the dose over 1 second to make up for clearance and 
        elimination loses starting at a point in time over a duration

        Parameters
        ----------
        infusion_list
            2d array of infusions with each row containing:
            [start time of infusion in seconds,
            dose of infusion over 1 second,
            duration of infusion in seconds,
            end time of infusion in seconds]
        target
            plasma/ effect site concentration to maintain at
        time
            start time in seconds to calculate maintenance infusion from
        duration
            time in seconds of the duration of the maintenance
            infusion

        Returns
        -------
        double
            dose per 1 second that would maintain a steady state
            concentration over the duration time period
        """

        cdef double dose, cp

        time += duration
        cp = self.calculate_cp(infusion_list, time)

        if target - cp <= 0:
            return 0

        if duration == 0:
            return 0  # avoid divide by zero error

        dose = (target - cp) / (
            (self.A / self.alpha) * (1 - exp(-self.alpha * duration)) +
            (self.B / self.beta) * (1 - exp(-self.beta * duration)) +
            (self.C / self.gamma) * (1 - exp(-self.gamma * duration)))

        return dose


    cpdef maintenance_infusion_list(self, double [:] target_concentration,
                                    double [:, :] infusion_list, int duration,
                                    int multiplier, double drug_concentration,
                                    int max_infusion_rate):

        """Returns the maintenance infusions required to maintain a target
        plasma/ effect site concentration over a time period

        Parameters
        ----------
        target_concentration
            1d array of targets with each row containing:
            [start time of target in seconds,
            target concentration,
            end time of target in seconds]
        infusion_list
            2d array of infusions with each row containing:
            [start time of infusion in seconds,
            dose of infusion over 1 second,
            duration of infusion in seconds,
            end time of infusion in seconds]
        duration
            time in seconds of the duration of each maintenance
            infusion
        multiplier
            multiplied to the duration on each iteration increasing
            the duration time

        Returns
        -------
        np.ndarray
            2d array of infusions with each row containing:
            [start time of infusion in seconds,
            dose of infusion over 1 second,
            duration of infusion in seconds,
            end time of infusion in seconds]
        """

        cdef int t, end_v

        inf_out = infusion_list
        start = int(target_concentration[0])
        target = target_concentration[1]
        end = int(target_concentration[2])

        # Initial infusion
        end_v = start + duration
        if end_v > end:
            end_v = end
            duration = end - start

        dose = self.maintenance_infusion(inf_out, target, start, duration)

        inf_v = np.array([start, dose, duration, end_v], dtype=np.float64)
        inf_out = np.vstack((inf_out, inf_v))

        # Remaining infusions
        t = start + duration
        duration *= multiplier

        while t < end:

            end_v = t + duration
            if end_v > end:
                end_v = end
                duration = end_v - t

            dose = self.maintenance_infusion(inf_out, target, t, duration)

            rate = (dose / drug_concentration) * (60 * 60)

            # If rate above max rate match max infusion rate
            if rate > max_infusion_rate:
                dose = rate / (60 * 60) * drug_concentration

            inf_v = np.array([t, dose, duration, end_v], dtype=np.float64)
            inf_out = np.vstack((inf_out, inf_v))

            t += duration
            duration *= multiplier
        return inf_out


    cpdef int plasma_decrement_time(self, int time, double target,
                                    double [:, :] infusion_list):

        """ Method returns time in seconds to reach a plasma target after
        all infusions are stopped.

        Parameters
        ----------
        time
            start time in seconds to calcuate decrement from
        target
            plasma site concentration to decrement to
        infusion_list
            2d array of infusions with each row containing:
            [start time of infusion in seconds,
            dose of infusion over 1 second,
            duration of infusion in seconds,
            end time of infusion in seconds]

        Returns
        -------
        int
            time in seconds to reach a plasma site target once all
            infusions are stopped
        """

        cdef int end, t, decrement_time
        cdef double cp
        cdef Py_ssize_t x_max = int(infusion_list.shape[0])
        cdef Py_ssize_t x

        inf_tmp = np.empty((0, 4), dtype=np.float64)

        for x in range(x_max):
            end =  int(infusion_list[x, 3])

            if end <= time:
                inf_tmp = np.vstack((inf_tmp, infusion_list[x, :]))
            else:
                inf_tmp = np.vstack((inf_tmp, infusion_list[x, :]))
                inf_tmp[x, 3] = time
                break

        t = time
        cp = self.calculate_cp(inf_tmp, time)

        if target == 0:
            target = 0.1  # Approximate 0 to 0.1 to avoid infinite loop

        while cp >= target:
            cp = self.calculate_cp(inf_tmp, t)
            t += 1

        decrement_time = t - time

        return decrement_time


    cpdef int effect_decrement_time(self, int time, double target,
                                    double [:, :] infusion_list):

        """ Method returns time in seconds to reach an effect target after
        all infusions are stopped.

        Parameters
        ----------
        time
            start time in seconds to calcuate decrement from
        target
            effect site concentration to decrement to
        infusion_list
            2d array of infusions with each row containing:
            [start time of infusion in seconds,
            dose of infusion over 1 second,
            duration of infusion in seconds,
            end time of infusion in seconds]

        Returns
        -------
        int
            time in seconds to reach a effect site target once all
            infusions are stopped
        """

        cdef int end, t, decrement_time
        cdef double previous_cp, current_cp, previous_ce
        cdef Py_ssize_t x_max = int(infusion_list.shape[0])
        cdef Py_ssize_t x

        previous_ce = 0
        inf_tmp = np.empty((0, 4), dtype=np.float64)

        for x in range(x_max):
            end =  int(infusion_list[x, 3])

            if end <= time:
                inf_tmp = np.vstack((inf_tmp, infusion_list[x, :]))
            else:
                inf_tmp = np.vstack((inf_tmp, infusion_list[x, :]))
                inf_tmp[x, 3] = time
                break

        t = 0
        current_cp = 0
        current_ce = 0

        if target == 0:
            target = 0.1  # Approximate 0 to 0.1 to avoid infinite loop

        while True:
            previous_cp = current_cp
            current_cp = self.calculate_cp(inf_tmp, t)

            current_ce = self.calculate_ce(current_cp, previous_cp,
                                           previous_ce)

            previous_ce = current_ce

            if (current_ce <= target) and (t > time):
                break

            t += 1

        decrement_time = t - time

        return decrement_time


    cpdef double ke0_tpeak_method_minimise(self, double ke0, double dose,
                                           double tpeak, double ce_tpeak):
        """ Minimisation function to solve ke0 'tpeak' method equation

        Parameters
        ----------
        ke0
            ke0 effect compartment equilibrium rate constant
        dose
            total dose of drug
        tpeak
            time in seconds of peak effect
        ce_tpeak
            effect site concentration at tpeak

        Returns
        -------
        double
            minimisation target

        """
        cdef double f = 0

        f = ((ke0 * self.A) / (ke0 - self.alpha)) * \
            (self.alpha * exp(-self.alpha * tpeak) - ke0 * exp(-ke0 * tpeak))

        f += ((ke0 * self.B) / (ke0 - self.beta)) * \
             (self.beta * exp(-self.beta * tpeak) - ke0 * exp(-ke0 * tpeak))

        f += ((ke0 * self.C) / (ke0 - self.gamma)) * \
             (self.gamma * exp(-self.gamma * tpeak) - ke0 * exp(-ke0 * tpeak))
        
        f *= dose
        f /= ce_tpeak

        return f


    cpdef double ke0_tpeak_method(self, double dose, double tpeak,
                                  double ce_tpeak):
        """ Returns ke0 using the 'tpeak' method equation

        Parameters
        ----------
        dose
            total bolus dose of drug given at time 0 over 1 second
        tpeak
            time in seconds of peak effect
        ce_tpeak
            effect site concentration at tpeak

        Returns
        -------
        double
            ke0 effect compartment equilibrium rate constant, units /second

        """
        root = optimize.brentq(self.ke0_tpeak_method_minimise,
                               a=1e-5, b=1e2,
                               args=(dose, tpeak, ce_tpeak))
        return root
