from src.hex.single_flow.flow_abc import HexFlow, FlowInputData, FlowOutletThermalData
from src.isobar.isobar import Isobar
from src.isobar.my_data_classes import Pressure, Temperature
# from flow_abc import FlowOutletThermalData

import numpy as np
import sys


class HexFlowHot(HexFlow):

    def __init__(self, inp_dat: FlowInputData, t_out_lim: Temperature):
        super().__init__(inp_dat, t_out_lim)

    def calc_out_thsdsdt_with_eps_qact(self, eps: float, q_act_w: float) -> FlowOutletThermalData:
        # q_act_w - общее кол-ва тепла, передаваемое от всех гор. потоков всем хол.потокам
        # eps - доля тепла данного потока в общей тепл. нагрузке: eps = flow[i].q_w / q_act_w
        q_w = -q_act_w * eps                                     # q_act_w > 0, но для гор. потока q_w < 0
        q_jmol = q_w / self.m_mols                               # q_w < 0 и q_jmol < 0
        h_out_jmol = self.inlet.h_jmol + q_jmol                       # h_out_jmol < self.inlet.h для гор. потока
        t_out_k, s_out_jmolk = self.calc_t_k_s_jmolk_with_h_via_isobar(h_jmol=h_out_jmol)
        dt_out_lim = t_out_k - self.outlet_lim.t_k                 # t_out_k > self.outlet_lim.t; hot flows only!!!

        # t_out_k control: self.outlet_lim.t_k задавалась, исходя из анализа тем-р встречных потоков на входе.
        # здесь t_out_k = f(h_out) и даже в том случае, когда h_out=self.outlet_lim.h_jmol возникает ошибка числ.
        # расчета t_out_k = f(h_out) и t_out_k может незначительно отличаться от self.outlet_lim.t_k на десятые градуса
        # ниже я пытаюсь искусственно исправить эту ситуацию
        if dt_out_lim < 0:
            if dt_out_lim < -0.5:
                print('t_out ', t_out_k, 't_out_lim ', self.outlet_lim.t_k, 'dt ', dt_out_lim)
                sys.exit('hot flow: t_out < t_out_lim in "calc_out_thsdsdt_with_eps_qact". crit.error')
            else:
                t_out_k = self.outlet_lim.t_k

        ds_jmolk = s_out_jmolk - self.inlet.s_jmolk                    # ds_jmolk < 0; hot flows only!!!
        ds_wk = ds_jmolk * self.m_mols

        return FlowOutletThermalData(t_k=t_out_k, h_jmol=h_out_jmol, s_jmolk=s_out_jmolk, dt=dt_out_lim, ds_wk=ds_wk)

    def calc_out_thsds_limiting_with_eps_qact(self, eps: float, eps_hex: float, q_act_w: float) -> [FlowOutletThermalData, float, float]:
        #
        # эта ф-ция "избыточна", т.е. она нужна лишь для случая, когда пакет потоков имеет: has_min_q_w=True
        # тогда каждый из пакетов потока имеет лишь по одному единств. значению eps и последующий расчет выход.
        # параметров потока, как f(q_act_w, eps) приводит к очень незначительному, но превышению выход. параметров
        # над их предел. значениями (просто из-за ошибки расчета св-в). чтобы этого избежать, после расчета eps
        # я пересчит. и заменяю имеющиеся предельные выходные параметры новыми, рассчит. как f(q_act_w, eps)
        # это происходит в ф-ции set_flow_eps
        # очевидно, что в этом случае нет нужды в сравнении выход. параметров с их предельными значениями, как это
        # имеет место в calc_out_thsdsdt_with_eps_qact
        #
        q_lim_w = -q_act_w/eps_hex * eps    # q_lim_w < 0 for hot flow, but q_act_w > 0;
                                            # q_act_w/eps_hex - возврат к изначальному значению sum(flow[i].q_w), кот.
                                            # могло быть больше q_act_w, т.к. q_act_w = sum(flow[i].q_w)*eps_hex
        q_lim_jmol = q_lim_w/self.m_mols
        h_out_lim_jmol = self.inlet.h_jmol + q_lim_jmol                      # "+" since q_lim_jmol < 0 for hot flow
        t_out_lim_k, s_out_lim_jmolk = self.calc_t_k_s_jmolk_with_h_via_isobar(h_jmol=h_out_lim_jmol)
        dt_out_lim = 0
        ds_lim_wk = (s_out_lim_jmolk - self.inlet.s_jmolk)*self.m_mols     # < 0; hot flows only!!!, [W/K]
        return FlowOutletThermalData(t_k=t_out_lim_k,
                                     h_jmol=h_out_lim_jmol,
                                     s_jmolk=s_out_lim_jmolk,
                                     dt=dt_out_lim,
                                     ds_wk=ds_lim_wk), q_lim_jmol, q_lim_w

    def calc_t_h_s_ds_distribution(self, t_out: float, h_out: float, s_out: float) -> [[], [], [], [], []]:
        # из super возвращаем t, h, s, dq в порядке: inlet - 0  -> outlet - self.n_nodes-1
        # ds в super не рассчитывался и возвращ. пустой список []
        t, h, s, ds, dq, dq_sum = super().calc_t_h_s_ds_distribution(t_out, h_out, s_out)  # return t, h, s, ds = [] без значений

        # условимся вести индексирование для всех потоков от входа хол. потока: inlet - 0  -> outlet - self.n_nodes-1
        # т.о. t, h, dq для  гор. потоков надо будет реверсировать:  outlet - 0  <- inlet - self.n_nodes-1
        ds = [0.0]
        for i in range(1, self.n_nodes):
            ds.append((s[i] - s[i - 1]) * self.m_mols)

        _t = list(reversed(t))
        _h = list(reversed(h))
        _s = list(reversed(s))
        _dq = list(reversed(dq))
        _ds = list(reversed(ds))
        dq_sum = _dq[0]
        return np.array(_t), np.array(_h), np.array(_s), np.array(_ds), np.array(_dq), dq_sum

    def section_flow_for_qt_calc(self, n_sections: int, n_nodes: int) -> None:
        super().section_flow_for_qt_calc(n_sections, n_nodes)

    def remove_needless_eps(self, _min: float, _max: float) -> None:
        super().remove_needless_eps(_min, _max)

    def set_eps_for_flows_in_active_pack(self, n_flows_in_pack: int, dlt_eps: float, q_act_w: float) -> []:
        eps_min = dlt_eps
        eps_max_basic = 1.0 - eps_min * (n_flows_in_pack - 1)  # calc. with n and dlt_eps ONLY
        n_epss = int(round((eps_max_basic - eps_min) / dlt_eps, 0)) + 1

        eps_list = []
        eps = eps_min
        for i in range(n_epss):
            q_w = q_act_w * eps
            q_jmol = q_w/self.m_mols
            h_out_jmol = self.inlet.h_jmol - q_jmol         # hot flow h_in > h_out
            # if self.isobar.h_min_jmol <= h_out_jmol < self.isobar.h_max_jmol:
            if h_out_jmol > self.outlet_lim.h_jmol:
                eps_list.append(round(eps, 4))
                eps += dlt_eps
            else:
                break
        return eps_list

    def set_eps_for_flows_in_passive_pack(self, eps_hex: float, q_act_w: float) -> []:
        return super().set_eps_for_flows_in_passive_pack(eps_hex, q_act_w)

    def set_isobar(self, t_in_k: float = None, t_out_k: float = None, p_kpa: float = None) -> Isobar:
        t_max_k = t_in_k + 10.0
        t_min_k = t_out_k - 10.0
        return Isobar(fluid=self.fluid,
                      p=Pressure(value=p_kpa, units='kpa'),
                      t_min=Temperature(value=t_min_k, units='k'),
                      t_max=Temperature(value=t_max_k, units='k'))

    def _set_mass_rate(self, m_value: float, m_units: str) -> (float, float):
        return super()._set_mass_rate(m_value, m_units)

    def _set_temperature(self, t_value: float, t_units: str) -> float:
        return super()._set_temperature(t_value, t_units)

    def _set_pressure(self, p_value: float, p_units: str) -> float:
        return super()._set_pressure(p_value, p_units)

    # calculate properties via isobar-----------------------------------------------------------------------------------
    #
    def calc_t_with_ph_via_isobar(self, h_jmol: float) -> float:
        return super().calc_t_with_ph_via_isobar(h_jmol)

    def calc_h_jmol_with_pt_via_isobar(self, t_k: float) -> float:
        return super().calc_h_jmol_with_pt_via_isobar(t_k)

    def calc_s_jmolk_with_pt_via_isobar(self, t_k: float) -> float:
        return super().calc_s_jmolk_with_pt_via_isobar(t_k)

    def calc_h_jmol_s_jmolk_with_pt_via_isobar(self, t_k: float) -> [float]:
        return super().calc_h_jmol_s_jmolk_with_pt_via_isobar(t_k)

    def calc_t_k_s_jmolk_with_h_via_isobar(self, h_jmol: float) -> [float]:
        return super().calc_t_k_s_jmolk_with_h_via_isobar(h_jmol=h_jmol)

    # calculate properties via rp10-------------------------------------------------------------------------------------
    #
    def calc_h_from_pt_via_rp10(self, p_kpa: float, t_k: float) -> (float, float):
        return super().calc_h_from_pt_via_rp10(p_kpa, t_k)

    def calc_t_from_ph_via_rp10(self, p_kpa: float, h_jmol: float) -> (float, float):
        return super().calc_t_from_ph_via_rp10(p_kpa, h_jmol)

    # def reverse_flow_cs_lists(self) -> None:
    #     pass
    #
    # def calc_dh_from_q(self, _q_w: float) -> float:
    #     # q = dh * n_s * m
    #     return super().calc_dh_from_q(_q_w)     # dh, [J/mol]; > 0 for both flows
    #
    # def get_t_as_list(self) -> list[float]:
    #     return super().get_t_as_list()
    #
    # def get_s_as_list(self) -> list[float]:
    #     return super().get_s_as_list()
