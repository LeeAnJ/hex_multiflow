from dataclasses import dataclass
from abc import ABC, abstractmethod
import sys
# import copy

from src.rp10.fluid.fluid_class import (RP10Fluid)
from src.isobar.isobar import Isobar
from src.isobar.my_data_classes import Pressure, Temperature, MassRate
from src.rp10.units_converters.units import UnitsInternal as Ui
from src.rp10.units_converters.units import UnitsUser as Uu
import src.rp10.units_converters.converters as conv


def set_isobar_for_hex_flow(fluid: RP10Fluid, p: Pressure, t_min: Temperature, t_max: Temperature) -> Isobar:
    # изобара должна работать в диапазоне [t_min .. t_max]. но, линейная интерполяция в numpy не работает в краевых
    # точках. для этого я ввожу dt_end = 10.0. Тогда температурный диапазон [(t_min - dt_end) .. (t_max + dt_end)]  и
    # в точках t_min, t_max линейная интерполяция будет работать нормально
    dt_end = 10.0
    _t_min = t_min.value - dt_end
    _t_max = t_max.value + dt_end
    return Isobar(fluid=fluid,
                  p=Pressure(value=p.value, units=p.units),
                  t_min=Temperature(value=_t_min, units=t_min.units),
                  t_max=Temperature(value=_t_max, units=t_max.units))


@dataclass
class FlowInputData:
    # fluid: RP10Fluid | None
    isobar: Isobar | None
    mass_rate: MassRate = MassRate(value=None, units=Uu.m)  # units user (Uu) for mass rate (m) Uu.m :[kgs]
    t_in: Temperature = Temperature(value=None, units=Ui.t)  # units internal (Ui) for tem-re (t) Ui.t :[k]
    p_in: Pressure = Pressure(value=None, units=Ui.p)  # units internal (Ui) for pressure (p) Ui.p :[kpa]

kelvin = float  # [K]
kpa = float  # [kPa]
jmol = float  # [J/mol]
jmolk = float  # [J/mol.K]
w = float  # [W]
wk = float  # [W/K]


@dataclass
class FlowThermalData:
    t_k: kelvin | None = None
    p_kpa: kpa | None = None
    h_jmol: jmol | None = None
    s_jmolk: jmolk | None = None


@dataclass(frozen=True)
class FlowThermalDataFrozen:
    t_k: kelvin | None = None
    p_kpa: kpa | None = None
    h_jmol: jmol | None = None
    s_jmolk: jmolk | None = None

# нельзя наследовать от FlowThermalData, т.к. FlowThermalData "не заморожен"
# нужно было создать такой же, но "замороженный" класс данных
# @dataclass(frozen=True)
# class FlowInletThermalData(FlowThermalDataFrozen):


@dataclass(frozen=True)
class FlowOutletThermalData(FlowThermalDataFrozen):
    dt: kelvin | None = None  # t_out_k - t_out_limiting_value: hot flows; t_out_limiting_value - t_out_k: cold flows
    ds_wk: wk | None = None  # ds = (s_out - s_in)*m_mols, j/mol.K * mol/s = J/s.K = W/K


@dataclass
class FlowOutletLimitingThermalData(FlowThermalData):
    dt: kelvin | None = None  # t_out_k - t_out_limiting_value: hot flows; t_out_limiting_value - t_out_k: cold flows
    ds_wk: wk | None = None  # ds = (s_out - s_in)*m_mols, j/mol.K * mol/s = J/s.K = W/K


class HexFlow(ABC):
    """
    self.m_mols, self.m_kgs
    self.t_in,            self.h_in,            self.p
    self.t_out_lim_value, self.h_out_lim_value, self.p
    self.q_jmol, self.q_w
    """

    @abstractmethod
    def __init__(self, inp_dat: FlowInputData, t_out_lim: Temperature):
        t_out_lim_k = conv.convert_arg_to_internal_units(t_out_lim.value, t_out_lim.units)  # t,[c] -> [k]
        # working fluid
        # self.fluid = inp_dat.fluid
        self.fluid = inp_dat.isobar.fluid
        self.isobar = inp_dat.isobar

        # mass rate in [mol/sec] and [kg/sec]:
        self.m_mols, self.m_kgs = self._set_mass_rate(inp_dat.mass_rate.value, inp_dat.mass_rate.units)

        # set t_in, p values in internal units (k, kpa) and isobar
        _t_in_k = self._set_temperature(inp_dat.t_in.value, inp_dat.t_in.units)
        self.p = self._set_pressure(inp_dat.p_in.value, inp_dat.p_in.units)  # pressure, [kPa]
        # self.isobar = self.set_isobar(t_in_k=_t_in_k, t_out_k=t_out_lim_k, p_kpa=self.p)

        # flow inlet port actual and invariable data: t, p, h, s will be frozen
        _h_in_jmol, _s_in_jmolk = self.isobar.get_h_jmol_s_jmolk_with_t__linear_interpolation(t_k=_t_in_k)
        self.inlet = FlowThermalDataFrozen(t_k=_t_in_k,
                                           p_kpa=self.p,
                                           h_jmol=_h_in_jmol,
                                           s_jmolk=_s_in_jmolk)

        # print(f' m_mols = {self.m_mols}, m_kgs = {self.m_kgs}')
        # print(f" inlet: t = {self.inlet.t_k} p = {self.inlet.p_kpa} h = {self.inlet.h_jmol} s = {self.inlet.s_jmolk}")

        # flow outlet port LIMITING data: t, p, h, s, ds, dt will be recalculated and are not frozen
        # here t_out_lim_k is __init__ argument calculated in hex_flows_init
        _h_out_lim_jmol, _s_out_lim_jmolk = self.isobar.get_h_jmol_s_jmolk_with_t__linear_interpolation(t_k=t_out_lim_k)
        self.outlet_lim = FlowOutletLimitingThermalData(t_k=t_out_lim_k,
                                                        p_kpa=self.inlet.p_kpa,
                                                        h_jmol=_h_out_lim_jmol,
                                                        s_jmolk=_s_out_lim_jmolk,
                                                        ds_wk=(_s_out_lim_jmolk - self.inlet.s_jmolk) * self.m_mols,
                                                        dt=0.0)

        # flow outlet port actual data are f(eps[i], q_act_w): to be calculated later with (eps[i], q_act_w) available
        self.outlet = FlowOutletThermalData()  # frozen t, p, h, s, ds, dt

        # для каждого еps из self.eps (см.ниже) проводим расчет self.outlet = FlowOutletThermalData()
        # и данные сохр. в списке self.outlet_dat_with_eps; потом по eps для потока находим t,h,s данные на выходе
        self.outlet_dat_with_eps: list[FlowThermalData] = list()

        # эти два параметра используются в self.set_eps_reduced() и в ф-ции пакета pack._calc_flows_total_q_w()
        self.q_lim_jmol = self.outlet_lim.h_jmol - self.inlet.h_jmol  # hot flows: < 0; cold flows > 0
        self.q_lim_w = self.q_lim_jmol * self.m_mols  # hot flows: < 0; cold flows > 0

        self.eps = []  # этот список теоретически здесь не нужен, он будет задан в pack.set_eps_per_flows ->
        # flow.set_eps_reduced или flow.set_eps_full_range, но к нему обращается ф-ция remove_needless_eps
        # так что я его оставляю

        # lists for qt-diagram data
        # будем начинать индексацию/нумерацию потоков с холодной стороны ТО:
        #  t_hot_out <- .......... <- t_hot_in
        #  0, 1, 2, ................. (n_nodes-1)
        #  t_cold_in -> .......... -> t_cold_out
        self.n_sects = None  # кол-во "секций", на которые разделен ТО, а значит и все потоки
        self.n_nodes = None  # кол-во "узлов" на секциях = кол-ву секций + 1
        self.node_index = []  # [0, 1, 2, ..., n_nodes-1]
        self.sect_index = []  # [0, 1, 2, ..., n_sects-1]

        # эти списки будут созданы и преобраз. в np.array в ф-ции calc_t_h_s_ds_distribution
        # self.t_k = []           # тем-ра потока в узлах секций: index E [0..self.n_nodes-1]
        # self.h_jmol = []        # энтальпия потока в узлах секций: index E [0..self.n_nodes-1]
        # self.s_jmolk = []       # энтропия потока в узлах секций: index E [0..self.n_nodes-1]

        # self.dh_jmol = []
        # self.dh_w = []
        # self.q_jmol = []   q[i] = h[i] - h[0], т.е. кол-во тепла наростает
        # self.q_w = []

        # self.n_cs = None  # number of cross-sections = number of elements in p,t,h lists
        # self.n_s = None  # number of sections = self.n_cs - 1
        # self.i_max = None  # maximal index in p,t,h lists (i_max = n_cs - 1)
        # self.cs_indexes: list[int] = list()  # [0, 1, 2, ... , self.i_ma] - list of indexes
        # self.cs: list[FlowThermalData] = list()  # list of p,t,h,s,ds_local,ds_integral for every cross-section
        # self.dh = None  # calc. in self.calc_dh_from_q(_q_w)

    @abstractmethod
    def calc_out_thsdsdt_with_eps_qact(self, eps: float, q_act_w: float) -> [FlowOutletThermalData, float, float]:
        pass

    @abstractmethod
    def calc_out_thsds_limiting_with_eps_qact(self,
                                              eps: float,
                                              eps_hex: float,
                                              q_act_w: float) -> [FlowOutletThermalData, float, float]:
        pass

    @abstractmethod
    def calc_t_h_s_ds_distribution(self, t_out: float, h_out: float, s_out: float) -> [[], [], [], [], []]:
        # будем вести индексирование от входа хол. потока. т.о. t, h для  гор. потока надо будет реверсировать !!!!
        # для хол. потока:  t[0] = t_in, h[0] = h_in  ... -> ... t[n] = t_out, h[n] = h_out
        # для гор. потока:  t[0] = t_out, h[0] = h_out  ... <- ... t[n] = t_in, h[n] = h_in
        t_in = self.inlet.t_k
        h_in = self.inlet.h_jmol
        s_in = self.inlet.s_jmolk
        dh = (h_out - h_in) / self.n_sects  # hot flow < 0; cold flow > 0
        h = [h_in + dh * i for i in self.node_index]
        # с нарастающим эффектом (т.е. dlt_q_w[self.node_index]-это все перед.тепло
        dq_w = [dh * self.m_mols * i for i in self.node_index]
        dq_w_sum = None  # пол-ное кол-во тепла будет в ф-циях hot_ cold_ потоков
        t = [t_in]
        s = [s_in]
        ds = []
        for i in range(1, self.n_sects):  # self.n_sects = self.n_nodes - 1, т.е. t_out надо добавить потом
            _t, _s = self.calc_t_k_s_jmolk_with_h_via_isobar(h_jmol=h[i])
            t.append(_t)
            s.append(_s)
        t.append(t_out)
        s.append(s_out)
        return t, h, s, ds, dq_w, dq_w_sum  # in flow_hot, flow_cold: calc ds, convert to np.array

    @abstractmethod
    def section_flow_for_qt_calc(self, n_sections: int, n_nodes: int) -> None:
        # ТО сексционируется на n_sections в зависимости от заданного dt_as (как правило dt_as = 10.0 K) в Hex __init__
        # после того, как иниц. гор.и хол. пакеты и установлены предельные гор. и хол. тем-ры ТО
        # задаем self.n_sects и self.n_nodes (n_nodes = n_sections + 1) после секц. собственно ТО
        self.n_sects = n_sections
        self.n_nodes = n_nodes
        self.node_index = [i for i in range(self.n_nodes)]
        self.sect_index = [i for i in range(self.n_sects)]

    @abstractmethod
    def remove_needless_eps(self, _min: float, _max: float) -> None:
        for eps in self.eps.copy():
            if eps < _min or eps > _max:
                self.eps.remove(eps)
        # зафиксировали "обрезанный" self.eps как кортеж. теперь его эл-ты изменять нельзя
        self.eps = tuple(self.eps)

    @abstractmethod
    def set_eps_for_flows_in_passive_pack(self, eps_hex: float, q_act_w: float) -> []:
        # этот вариант для пакета с меньшим суммарным теплосодержанием потоков.
        # для него мне известно предельное теплосодержание каждого из потоков: q_lim = f(t_out=t_out_lim)
        # в дальнейшем q_act = sum(q_lim)*eps_hex. т.о. я смогу передать в этот пакет чуть меньше тепла, чем
        # он теоретически мог бы принять. при этом я не знаю, как это тепло распределится между потоками пакета.
        # пока буду считать, что сохранится начальное распределение: eps[i] = q_lim[i]/sum(q_lim). т.е. уменьшение
        # sum(q_lim) на eps_hex привело также к пропорц. уменьшению каждого из q_lim[i], а отношение осталось неизменным

        _sum_q_lim_w = q_act_w / eps_hex  # возврат от q_act_w к sum(q_lim), т.к. в экз.потока мне неизвест. sum(q_lim)
        eps = tuple([round(abs(self.q_lim_w) / _sum_q_lim_w, 4), ])  # for cold flow self.q_lim_w < 0
        return eps

    @abstractmethod
    def set_eps_for_flows_in_active_pack(self, n_flows_in_pack: int, dlt_eps: float, q_act_w: float) -> []:
        pass

    @abstractmethod
    def set_isobar(self, t_in_k: float = None, t_out_k: float = None, p_kpa: float = None) -> Isobar:
        pass

    @abstractmethod
    def _set_mass_rate(self, m_value: float, m_units: str) -> (float, float):
        _m_units = conv.convert_units_string(m_units)
        if _m_units == Ui.m.value:  # mass rate input in internal units
            m_mols = m_value
            m_kgs = conv.convert_mass_rate_to_user_units(m_value,
                                                         _m_units,
                                                         self.fluid.mm_g_mol)
        elif _m_units == Uu.m.value:  # mass rate input in user units
            m_kgs = m_value
            m_mols = conv.convert_mass_rate_to_internal_units(m_value,
                                                              _m_units,
                                                              self.fluid.mm_g_mol)
        else:
            sys.exit('mass rate units neither internal nor user in "set_mass_rate". Program terminated')
        return m_mols, m_kgs

    @abstractmethod
    def _set_temperature(self, t_value: float, t_units: str) -> float:
        _t_units = conv.convert_units_string(t_units)
        if _t_units == Ui.t.value:  # temperature input in internal units
            t = t_value
        elif _t_units == Uu.t.value:  # temperature rate input in user units
            t = conv.convert_arg_to_internal_units(t_value, _t_units)
        else:
            sys.exit('temperature units neither internal nor user in "set_temperature". Program terminated')
        return t

    @abstractmethod
    def _set_pressure(self, p_value: float, p_units: str) -> float:
        _p_units = conv.convert_units_string(p_units)
        if _p_units == Ui.p.value:  # pressure input in internal units
            p = p_value
        elif _p_units == Uu.p.value:  # pressure rate input in user units
            p = conv.convert_arg_to_internal_units(p_value, _p_units)
        else:
            sys.exit('pressure units neither internal nor user in "set_pressure". Program terminated')
        return p

    # calculate properties via isobar-----------------------------------------------------------------------------------
    #
    @abstractmethod
    def calc_t_with_ph_via_isobar(self, h_jmol: float) -> float:
        return self.isobar.get_t_k_with_h_jmol__linear_interpolation(h_jmol=h_jmol)

    @abstractmethod
    def calc_h_jmol_with_pt_via_isobar(self, t_k: float) -> float:
        return self.isobar.get_h_jmol_with_t__linear_interpolation(t_k=t_k)

    @abstractmethod
    def calc_s_jmolk_with_pt_via_isobar(self, t_k: float) -> float:
        return self.isobar.get_s_jmolk_with_t__linear_interpolation(t_k=t_k)

    @abstractmethod
    def calc_h_jmol_s_jmolk_with_pt_via_isobar(self, t_k: float) -> [float]:
        return self.isobar.get_h_jmol_with_t__linear_interpolation(t_k=t_k), \
            self.isobar.get_s_jmolk_with_t__linear_interpolation(t_k=t_k)

    @abstractmethod
    def calc_t_k_s_jmolk_with_h_via_isobar(self, h_jmol: float) -> [float]:
        return self.isobar.get_t_k_s_jmolk_with_h__linear_interpolation(h_jmol=h_jmol)  # t, s = self.isobar.get_t_k_s_

    # calculate properties via rp10-------------------------------------------------------------------------------------
    #
    @abstractmethod
    def calc_h_from_pt_via_rp10(self, p_kpa: float, t_k: float) -> float:
        self.fluid.calc_spec_state(t=(t_k, 'k'), p=(p_kpa, 'kpa'))
        if self.fluid.error.index > 0:
            self.fluid.error.print_and_terminate()
        else:
            return self.fluid.state.get_data(flag='blk', x_symbol='h', x_units='jmol')

    @abstractmethod
    def calc_t_from_ph_via_rp10(self, p_kpa: float, h_jmol: float) -> (float, float):
        self.fluid.calc_spec_state(h=(h_jmol, 'jmol'), p=(p_kpa, 'kpa'))
        if self.fluid.error.index > 0:
            self.fluid.error.print_and_terminate()
        else:
            return self.fluid.state.get_data(flag='blk', x_symbol='t', x_units='k')

    # устаревшие ф-ции; пока не используются --------------------------------------------------------------------------
    #
    # @abstractmethod
    # def reverse_flow_cs_lists(self) -> None:
    #     for i in range(self.n_cs):
    #         last_item = self.cs.pop()
    #         self.cs.insert(i, last_item)
    #
    # @abstractmethod
    # def calc_dh_from_q(self, _q_w: float) -> float:
    #     # q = dh * n_s * m
    #     return _q_w / self.n_s / self.m_mols  # dh, [J/mol]; > 0 for both flows
    #
    # @abstractmethod
    # def get_t_as_list(self) -> list[float]:
    #     _t = list()
    #     for i in self.cs_indexes:
    #         _t.append(self.cs[i].t)
    #     return _t
    #
    # @abstractmethod
    # def get_s_as_list(self) -> list[float]:
    #     _s = list()
    #     for i in self.cs_indexes:
    #         _s.append(self.cs[i].s)
    #     return _s
