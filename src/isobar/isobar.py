import sys
import numpy as np

from src.isobar.my_data_classes import Pressure, Temperature, PthDat
from src.rp10.units_converters import converters as conv
from src.rp10.fluid.fluid_class import RP10Fluid


class Isobar:
    def __init__(self,
                 fluid: RP10Fluid,
                 p: Pressure,
                 t_min: Temperature,
                 t_max: Temperature) -> None:

        self.fluid = fluid
        self.p_kpa = conv.convert_arg_to_internal_units(p.value, p.units)  # p, [kPa]

        self.t_min_k = conv.convert_arg_to_internal_units(t_min.value, t_min.units)  # t, [k]
        self.t_max_k = conv.convert_arg_to_internal_units(t_max.value, t_max.units)  # t, [k]

        # calc. t,h on saturated line: bubble and dew points
        self.t_bubble_k, self.h_bubble_jmol, self.s_bubble_jmolk = self._bubble_point()
        self.t_dew_k, self.h_dew_jmol, self.s_dew_jmolk = self._dew_point()

        # init t, h and s arrays: just 2 points - [t_min, t_max], [h_min, h_max]  and [s_min, s_max]
        self.array_n, self.t_k, self.h_jmol, self.s_jmolk = self.init__t_h_s_arrays()

        self.h_min_jmol = self.h_jmol[0]  # used in recursion
        self.h_max_jmol = self.h_jmol[-1]  # used in recursion

        # at present entropy is not used as an argument for interpolation procedure
        # self.s_min_jmolk = self.s_jmolk[0]  # used in recursion
        # self.s_max_jmolk = self.s_jmolk[-1]  # used in recursion

        # insert t,h,s for bubble and dew points into self.t_k, self.h_jmol, self.s_jmol - arrays
        self.insert_ths_into_arrays(t_k=self.t_bubble_k, h_jmol=self.h_bubble_jmol, s_jmolk=self.s_bubble_jmolk)
        self.insert_ths_into_arrays(t_k=self.t_dew_k, h_jmol=self.h_dew_jmol, s_jmolk=self.s_dew_jmolk)

        # calc. additional t,h points in recurs. func and insert these points into: self.t_k, self.h_jmol
        self.calc_t_recursively()

    def calc_t_recursively(self) -> None:
        _inp = [[self.p_kpa, self.t_min_k, self.h_min_jmol]]

        if self.t_min_k < self.t_bubble_k:  # isobar starts in liq
            if self.t_bubble_k < self.t_max_k < self.t_dew_k:  # isobar ends in 2ph
                _inp.append([self.p_kpa, self.t_bubble_k, self.h_bubble_jmol])
            elif self.t_max_k > self.t_dew_k:  # isobar ends in vap
                _inp.append([self.p_kpa, self.t_bubble_k, self.h_bubble_jmol])
                _inp.append([self.p_kpa, self.t_dew_k, self.h_dew_jmol])
        elif self.t_bubble_k < self.t_min_k < self.t_dew_k:  # isobar starts in 2ph or vap
            if self.t_max_k > self.t_dew_k:  # isobar ends in vap
                _inp.append([self.p_kpa, self.t_dew_k, self.h_dew_jmol])
        _inp.append([self.p_kpa, self.t_max_k, self.h_max_jmol])

        for i in range(len(_inp) - 1):
            self.split_dt_recursively(inlet=PthDat(p=_inp[i][0], t=_inp[i][1], h=_inp[i][2]),
                                      outlet=PthDat(p=_inp[i+1][0], t=_inp[i+1][1], h=_inp[i+1][2]))

    # рекурс. функцию нельзя "кэшировать" (@functools.cache), т.к. аргументы типа PthDat - это не хэшируемый тим данных
    def split_dt_recursively(self, inlet: PthDat, outlet: PthDat) -> None:
        # linear interpol.:
        t_middle = (inlet.t + outlet.t) / 2.0
        h_middle = (outlet.h - inlet.h) / (outlet.t - inlet.t) * (t_middle - inlet.t) + inlet.h
        # rp10 calc.:
        h_middle_rp10, s_middle_rp10 = self.get_h_jmol_s_jmolk_with_t__rp10(t_k=t_middle)
        self.insert_ths_into_arrays(t_k=t_middle, h_jmol=h_middle_rp10, s_jmolk=s_middle_rp10)
        err = abs(h_middle_rp10 - h_middle) / h_middle_rp10 * 100.0

        if err > 0.5:
            self.split_dt_recursively(inlet=PthDat(p=self.p_kpa, t=inlet.t, h=inlet.h),
                                      outlet=PthDat(p=self.p_kpa, t=t_middle, h=h_middle_rp10))
            self.split_dt_recursively(inlet=PthDat(p=self.p_kpa, t=t_middle, h=h_middle_rp10),
                                      outlet=PthDat(p=self.p_kpa, t=outlet.t, h=outlet.h))
        else:
            return

    def init__t_h_s_arrays(self) -> (int, np.array, np.array, np.array):
        """
        create 2 arrays of 2 elements only: [t_min, t_max]; [h_min, h_max]
        :return: n-length of arrays; t_k, h_jmol
        """
        _t_k = np.array([self.t_min_k, self.t_max_k])
        _h_min_jmol, _s_min_jmolk = self.get_h_jmol_s_jmolk_with_t__rp10(t_k=self.t_min_k)
        _h_max_jmol, _s_max_jmolk = self.get_h_jmol_s_jmolk_with_t__rp10(t_k=self.t_max_k)
        _h_jmol = np.array([_h_min_jmol, _h_max_jmol])
        _s_jmolk = np.array([_s_min_jmolk, _s_max_jmolk])
        return len(_t_k), _t_k, _h_jmol, _s_jmolk

    def insert_ths_into_arrays(self, t_k: float, h_jmol: float = None, s_jmolk: float = None) -> None:
        if self.t_is_within_tmin_tmax(t_k):
            if h_jmol is None or s_jmolk is None:
                h_jmol, s_jmolk = self.get_h_jmol_s_jmolk_with_t__rp10(t_k=t_k)
            _index = np.searchsorted(self.t_k, t_k, side='right', sorter=None)
            self.array_n += 1
            self.t_k = np.insert(self.t_k, _index, t_k)
            self.h_jmol = np.insert(self.h_jmol, _index, h_jmol)
            self.s_jmolk = np.insert(self.s_jmolk, _index, s_jmolk)

    def t_is_within_tmin_tmax(self, t_k: float) -> bool:
        flag = False
        if self.t_min_k < t_k < self.t_max_k:
            flag = True
        return flag

    def _bubble_point(self) -> (float, float):  # t_bubble_k, h_bubble_jmol
        self.fluid.calc_sat_state(sat_curve_flag='l', p=(self.p_kpa, 'kpa'))

        if self.fluid.error.index > 0:
            self.fluid.error.print_and_terminate()
        _t_l_k = self.fluid.state.get_data(flag='bubble', x_symbol='t', x_units='k')
        _h_l_jmol = self.fluid.state.get_data(flag='bubble', x_symbol='h', x_units='jmol')
        _s_l_jmolk = self.fluid.state.get_data(flag='bubble', x_symbol='s', x_units='jmolk')
        return _t_l_k, _h_l_jmol, _s_l_jmolk

    def _dew_point(self) -> (float, float):  # t_dew_k, h_dew_jmol
        self.fluid.calc_sat_state(sat_curve_flag='v', p=(self.p_kpa, 'kpa'))

        if self.fluid.error.index > 0:
            self.fluid.error.print_and_terminate()
        _t_v_k = self.fluid.state.get_data(flag='dew', x_symbol='t', x_units='k')
        _h_v_jmol = self.fluid.state.get_data(flag='dew', x_symbol='h', x_units='jmol')
        _s_v_jmolk = self.fluid.state.get_data(flag='dew', x_symbol='s', x_units='jmolk')
        return _t_v_k, _h_v_jmol, _s_v_jmolk

    # h,s = f(t) along isobar.p_kpa
    def get_h_jmol_s_jmolk_with_t__linear_interpolation(self, t_k: float) -> [float]:
        if t_k < self.t_min_k or t_k > self.t_max_k:
            sys.exit('argument t_k in "get_h_jmol_s_jmolk_with_linear_interpolation" is out of acceptable range')
        i = np.searchsorted(self.t_k, t_k, side="right")
        dt = self.t_k[i]-self.t_k[i-1]
        dt_ = t_k-self.t_k[i-1]
        h = (self.h_jmol[i] - self.h_jmol[i-1])/dt * dt_ + self.h_jmol[i-1]
        s = (self.s_jmolk[i] - self.s_jmolk[i-1])/dt * dt_ + self.s_jmolk[i-1]
        return h, s

    # t,s = f(h) along isobar.p_kpa
    def get_t_k_s_jmolk_with_h__linear_interpolation(self, h_jmol: float) -> [float]:
        if h_jmol < self.h_min_jmol or h_jmol > self.h_max_jmol:
            sys.exit('argument h_jmol in "get_t_k_s_jmolk_with_h__linear_interpolation" is out of acceptable range')
        i = np.searchsorted(self.h_jmol, h_jmol, side="right")
        dh = self.h_jmol[i]-self.h_jmol[i-1]
        dh_ = h_jmol-self.h_jmol[i-1]
        t = (self.t_k[i] - self.t_k[i-1])/dh * dh_ + self.t_k[i-1]
        s = (self.s_jmolk[i] - self.s_jmolk[i-1])/dh * dh_ + self.s_jmolk[i-1]
        return t, s

    # s = f(t) along isobar.p_kpa
    def get_s_jmolk_with_t__linear_interpolation(self, t_k: float) -> float:
        if t_k < self.t_min_k or t_k > self.t_max_k:
            sys.exit('argument t_k in "get_s_jmolk_with_t__linear_interpolation" is out of acceptable range')
        return np.interp(t_k, self.t_k, self.s_jmolk)

    # h = f(t) along isobar.p_kpa
    def get_h_jmol_with_t__linear_interpolation(self, t_k: float) -> float:
        if t_k < self.t_min_k or t_k > self.t_max_k:
            sys.exit('argument t_k in "get_h_jmol_with_t__linear_interpolation" is out of acceptable range')
        return np.interp(t_k, self.t_k, self.h_jmol)

    # t = f(h) along isobar.p_kpa
    def get_t_k_with_h_jmol__linear_interpolation(self, h_jmol: float) -> float:
        if h_jmol < self.h_min_jmol or h_jmol > self.h_max_jmol:
            sys.exit('argument h_jmol in "get_t_k_with_h_jmol__linear_interpolation" is out of acceptable range')
        return np.interp(h_jmol, self.h_jmol, self.t_k)

# ---------------------------------- RefProp 10.0 ----------------------------------------------------------------------
    # h,s = f(t) with RefProp10
    def get_h_jmol_s_jmolk_with_t__rp10(self, t_k: float) -> [float]:
        self.fluid.calc_spec_state(t=(t_k, 'k'), p=(self.p_kpa, 'kpa'))
        if self.fluid.error.index > 0:
            self.fluid.error.print_and_terminate()
        else:
            return self.fluid.state.get_data(flag='blk', x_symbol='h', x_units='jmol'), \
                   self.fluid.state.get_data(flag='blk', x_symbol='s', x_units='jmolk')

    # t = f(t) with RefProp10
    def get_t_k_with_refprop10(self, h_jmol: float) -> float:
        self.fluid.calc_spec_state(h=(h_jmol, 'jmol'), p=(self.p_kpa, 'kpa'))
        if self.fluid.error.index > 0:
            self.fluid.error.print_and_terminate()
        else:
            return self.fluid.state.get_data(flag='blk', x_symbol='t', x_units='k')

    # t,s = f(h) with RefProp10
    def get_t_k_s_jmolk_with_h__rp10(self, h_jmol: float) -> [float]:
        self.fluid.calc_spec_state(h=(h_jmol, 'jmol'), p=(self.p_kpa, 'kpa'))
        if self.fluid.error.index > 0:
            self.fluid.error.print_and_terminate()
        else:
            return self.fluid.state.get_data(flag='blk', x_symbol='t', x_units='k'), \
                   self.fluid.state.get_data(flag='blk', x_symbol='s', x_units='jmolk')
