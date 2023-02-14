import sys
from dataclasses import dataclass, field
import math
import typing
import numpy as np
from itertools import product
from typing import Literal

from src.hex.single_flow.flow_cold import HexFlowCold
from src.hex.single_flow.flow_hot import HexFlowHot
from src.hex.errors import terminate_program
from src.hex.pack_of_flows.pack_hot import HexPackOfHotFlows
from src.hex.pack_of_flows.pack_cold import HexPackOfColdFlows

kelvin = float
percent = float
PackFlag = Literal["hot", "cold"]


@dataclass
class HexData:
    dt_as: kelvin  # max dt along hex section (i..i+1)
    dt_cs: kelvin  # min dt across section (between the hex's flows)
    # dt_error: kelvin  # acceptable range of dt_min E [dt_cs-dt_error .. dt_cs+dt_error]
    eps_hex: float
    dlt_eps: float
    # dlt_eps_hot_flow: float
    # dlt_eps_cold_flow: float
    # eps_error: percent  # acceptable range of eps: (1.0 - eps_min/eps_max)*100.0 <= eps_error


@dataclass
class HexFlows:
    hot: list[HexFlowHot]
    cold: list[HexFlowCold]


@dataclass
class HexPack:
    hot: HexPackOfHotFlows
    cold: HexPackOfColdFlows
    active: HexPackOfHotFlows | HexPackOfColdFlows | None
    passive: HexPackOfHotFlows | HexPackOfColdFlows | None


class Hex:

    def __init__(self, flows_sorted_with_tin: HexFlows, data: HexData):
        # eps_hex: характер. общие потери тепла при передаче от потоков гор. пакета ТО потокам хол.:
        # eps_hex = Qact/Qmax, где Qmax=min(sum(Qhot[i]), sum(Qcold[i]))
        self.eps_hex = data.eps_hex

        # организуем заданные холод. и гор. потоки в соответствующие пакеты: self.pack.hot, self.pack.cold
        # в пакетах потоки отсортированы по температурам на выходе (см. hot/cold_flows_inp)
        # self.pack.hot = HexPackOfHotFlows(flows=flows_sorted_with_tin.hot, dlt_eps=data.dlt_eps)
        # self.pack.cold = HexPackOfColdFlows(flows=flows_sorted_with_tin.cold, dlt_eps=data.dlt_eps)
        self.pack = HexPack(hot=HexPackOfHotFlows(flows=flows_sorted_with_tin.hot, dlt_eps=data.dlt_eps),
                            cold=HexPackOfColdFlows(flows=flows_sorted_with_tin.cold, dlt_eps=data.dlt_eps),
                            active=None, passive=None)

        # определим пакет с мин.теплосодержанием (он будет пассивным, т.е. с неизменяющимся набором eps для потоков)
        # предельное кол-ва тепла, передаваемое в ТО = мин.теплосодержанию=self.q_lim_w;
        # фактичесое кол-ва тепла, передаваемое в ТО (т.е. с учетом потерь): self.q_act_w = self.q_lim_w * self.eps_hex
        self.q_lim_w, self.q_act_w, \
            self.pack.active, self.pack.passive, \
            self.pack.hot.has_min_q_w, self.pack.cold.has_min_q_w = self._set_pack_with_lim_heat_capacity(
            eps_hex=self.eps_hex)

        # только после того как пакеты инициированы:
        # ТО по длине разделяется на self.n_qt_n_sect условных секций; при этом образуется self.n_qt_nods "узлов",
        # в кот. задаются-определяются t,h параметры потоков. кол-во узлов зависит от значений предельных тем-р в ТО и
        # "желаемого" шага по температуре по длине секции: self.dt_sect
        # после "секционирования" ТО задаем кол-ва секций и узлов для каждого потока в пакете (чтобы потом не возиться)
        self.dt_sect = data.dt_as
        self.dt_inter_flows_min = data.dt_cs  # min dt across section (between the hex's flows)
        self.n_qt_sect, self.n_qt_nods = self.split_into_sections_with_dt()  # self.dt_sect уже должно быть задано
        self._set_sections_in_flows(n_sections=self.n_qt_sect, n_nodes=self.n_qt_nods)

        """" 
        1. создадим СПИСКИ [eps,] для каждого потока в обоих пакетах: от eps_min=dlt_eps с шагом dlt_eps до 
           eps_max=f(n_flows, dlt_eps) при условии, что Qact*eps_max <= Qlim для данного потока. т.к. Qlim для потоков
           могут отличаться, то и eps_max для разных потоков (в активном пакете) могут ИМЕТЬ РАЗНЫЕ ЗНАЧЕНИЯ!!!
           списки [eps,] для потоков хранятся в:
                                self.pack.active.flow[i].eps, i E [1,n_flows]  на этом этапе это СПИСОК: list
                                self.pack.passive.flow[i].eps, i E [1,n_flows]  на этом этапе это уже КОРТЕЖ: tuple
            !!! если в пакете только один поток, то в списке [eps,] будет только одно значение eps=1: [1] !!!
            !!! для потоков активного пакета списки [eps,] практически всегда будут отличаться значением eps_max !!!
            !!! для потоков пассив.пакета в списках [eps,] будет только по 1му значению: [eps=Q[indx_потока]/Qact] !!!
         2. скомбинируем [eps,] потоков каждого пакета в tuples of tuples ((eps, eps),(eps, eps),...), так 
            чтобы sum(eps) каждого внутр. tuple = 1:
                                self.pack.active/passive.eps_sets.
            Смысл такого комбинирования в том, что каждый внутренний tuple: (eps, eps) является как бы отдельным
            "рабочим режимом" ТО, для которого можно рассчитать свои выходные параметры потоков, распределения
            температур и провести оценку работоспособности qt-диаграммы.
         3. т.к. eps_max в списках self.pack.active.flow[i].eps могут отличаться, то и при комбинировании eps c учетом
            равенства их суммы 1, некоторые минимальные значения eps из исходных списков в эти внутр. tuple не попали,
            и, следовательно, они больше не нужны. "обрезаем" лишние eps из списков eps для потоков активного пакета и 
            превращаем эти списки в tuples:
                                self.pack.active/passive.flow[i].eps, i E [1,n_flows]  на этом этапе это уже tuples
        """
        self._set_eps_combinations_useful_for_operating_modes()

        """" 
        1. расчет FlowOutletThermalData: frozen t, p, h, s, ds, dt данных в выходном сечении каждого потока для  
        каждого значения eps из кортежей self.pack.active/passive.flow[i]. рез-ты расчетов запишем в 2e матрицы-кортежи: 
        self.pack.active/passive.flows_outlet_dat[i_flow][j_eps]
        2. скомбинируем выходные данные для потоков каждого пакета в кортежи кортежей, аналогичные по структуре
        self.pack.active/passive.eps_sets; т.е. каждый внутренний кортеж будет аналогом некоторого отдельного 
        независимого"рабочего режима" ТО, в котором известны все выходные данные потоков:
        self.pack.active/passive.flows_outdat_in_sets[i_operating_mode][j_flow] -> j-flow's outlet t, p, h, s, ds, dt
        """
        self.pack.active.flows_outlet_dat, \
            self.pack.passive.flows_outlet_dat = self._calc_flows_outlet_dat_with_eps()

        self.pack.active.flows_outlet_dat_in_sets, \
            self.pack.passive.flows_outlet_dat_in_sets = self._set_flows_outdat_combinations_useful_for_operating_modes()

        # разница тем-р между самым хол. потоком в гор.пакете и самым гор. потоком в хол.пакете
        self.dt_inter_pack = []
        self.dt_inter_flows = None
        self.dt_min_inter_flows = None
        self.dt_max_inter_flows = None
        self.dt_avr_inter_flows = None
        self.ds_wk = []
        self.ds_wk_sum = None

        # for pack in [self.pack.active, self.pack.passive]:
        # print(pack.flows_outlet_dat)
        # print(pack.flows_outlet_dat_in_sets)
        # for i in range(pack.n_flows):
        #     print(i, pack.flows_outlet_dat[i])

        # self.ds_with_operating_modes_sorted = self.calc_ds_with_operating_modes()
        # print(self.ds_with_operating_modes_sorted)
        #
        # sys.exit()

        # print(self.pack.hot.flows_outlet_dat)
        # print(self.pack.cold.flows_outlet_dat)
        # print(self.pack.active.flows_outlet_dat)
        # print(self.pack.passive.flows_outlet_dat)
        # sys.exit('---')

        # выше:
        # 1. мы скомбинировали списки [eps] для каждого потока в обеих пакетах АКТИВНОГО пакета в кортеж кортежей таким образом, чтобы суммы eps в
        # мы скомбинировали списки [eps] для каждого потока в обеих пакетах АКТИВНОГО пакета в кортеж кортежей таким образом, чтобы суммы eps в
        # каждом кортеже = 1; для всех eps всех потоков в каждом пакете провели расчет выходных данных.
        # теперь соберем аналогичные по структуре кортежи кортежей (для каждого пакета), но уже из выходных данных.
        # т.е. в этом кортеже вместо eps[i][j] будет стоять соответств. out_dat[i][j].
        #
        # ВНИМАНИЕ:

        # self.pack.hot.flows_outlet_dat_in_sets = None
        # self.pack.cold.flows_outlet_dat_in_sets = None

        # здесь мы собрали выход. данные в tuple of tuples в соответствии с наборами pack.hot/c.eps_sets, полученными в
        # self._combine_flows_epss_in_tuple_of_eps_sets() и отсортировали их в зависимости от произв.энтропии
        # ВНИМАНИЕ! после сортировки наборы pack.hot/c.flows_outdat_in_sets могут уже не соответствовать наборам
        # pack.hot/c.eps_sets
        # self.set_flows_outdat_in_sets()

        # print(self.pack.hot.flows_outdat_in_sets)
        # print(self.pack.cold.flows_outdat_in_sets)
        # sys.exit('----')

    # ------------------------------------------------------------------------------------------------------------------
    # METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--METHODS--------
    # ------------------------------------------------------------------------------------------------------------------
    def calc_ds_with_operating_modes(self):
        ds_wk = []  # здесь будем хранить кортежи с инд.рабочего режима ТО и произведенной в этом реж. энтропии

        # the passive pack: пассивном пакете только 1 набор выходных данных по потокам, т.к. только один набор eps_set
        ds_sum_passive = 0.0
        for data in (self.pack.passive.flows_outlet_dat_in_sets[0]):
            ds_sum_passive += data.ds_wk

        # the active pack:
        for i_mode, data_set in enumerate(self.pack.active.flows_outlet_dat_in_sets):
            ds_sum = ds_sum_passive  # чтобы сразу суммировать энтроп. произвед в пассив. и актив. пакетах
            for i_flow, data in enumerate(data_set):
                ds_sum += data.ds_wk  # это уже значение энтропии, сгенерированной в ТО в данном рабочем режиме
            ds_wk.append(tuple([i_mode, ds_sum]))

        # отсортируем кортежи x = (i_mode, ds) по возрастанию энтропии (кот. записана во втором эл-те кортежа: x[1])
        ds_wk.sort(key=lambda x: x[1])
        # вернем как кортеж кортежей, т.к. изменять его больше не понадобится
        return tuple(ds_wk)

    def _set_eps_combinations_useful_for_operating_modes(self):
        # задаем для КАЖДОГО ИЗ ПОТОКОВ обоих пакетов массивы eps (в виде СПИСКОВ-list):
        #                         self.pack.passive/active.flow[i].eps, где i E [1..n_flows]
        # Для "пассивного" пакета в списке eps для каждого потока будет лишь одно значение. Для "активного" пакета
        # каждый поток обладает массивом значений eps [eps_min=dlt_eps .. eps_max=f(n_flows)].
        # Т.к. для потока должно выполняться условие: h_out_lim(eps_max, q_act, m_mols) не должно превышать предельного
        # значения (определяется ранее), то в активном пакете потоки имеют одинаковые eps_min=dlt_eps,
        # НО ВОЗМОЖНЫ РАЗНЫЕ EPS_MAX!!! Списки будут нуждаться в уточнении
        self.pack.passive.set_eps_per_flows(q_act_w=self.q_act_w, eps_hex=self.eps_hex)
        self.pack.active.set_eps_per_flows(q_act_w=self.q_act_w, eps_hex=self.eps_hex)

        # for pack in [self.pack.active, self.pack.passive]:
        #     # print(pack.eps_sets)
        #     for i in range(pack.n_flows):
        #         print(i, pack.flow[i].eps)
        # sys.exit()

        # собрать все возможные комбинации eps потоков данного пакета, которые дают в сумме 1.0, в список с эл-ми tuple,
        # где каждый элемент списка - tuple - это набор eps, дающий распределение q_act_w по потокам текущего пакета
        # для потоков "пассивного" пакета будет лишь одна комбинация, т.к. каждый поток имеет лишь одно значение eps
        # для потоков "активного" пакета таких комбинаций (tuple) будет много, но некоторые значения eps для отдельных
        # потоков пакета не войдут в эти комбинации, т.к. в сумме они не позволят получить 1. Т.о., эти "лишние" eps
        # надо в дальнейшем убрать из начальных списков eps для потоков, чтобы не мешались под ногами.
        self.pack.passive.eps_sets = tuple(self.pack.passive.combine_flows_epss_in_eps_sets())
        self.pack.active.eps_sets = tuple(self.pack.active.combine_flows_epss_in_eps_sets())

        # зная eps_min и eps_max для каждого потока, "обрежем" массивы pack.active.flow[i].eps для каждого потока
        # эти обрезанные массивы нужны для дальнейших расчетов выход.параметров потоков в зависимости от
        # текущих еps: h_out = f(q_act_w * eps[j]) -> t_out, s_out = f(h_out)
        self.pack.active.eps_min, self.pack.active.eps_max = self.pack.active.get_eps_min_max_from_pack_eps_combined()
        self.pack.active.from_flows_remove_needless_eps(_mins=self.pack.active.eps_min, _maxs=self.pack.active.eps_max)

    def _set_pack_with_lim_heat_capacity(self, eps_hex: float):
        # self.pack.hot.sum_of_q_lim_w < 0        гор. потоки "отдают" тепло
        # self.pack.cold.sum_of_q_lim_w > 0        хол. потоки "принимают" тепло

        pack_hot_has_min_q_w = False
        pack_cold_has_min_q_w = False

        q_lim_w = min([abs(self.pack.hot.sum_of_q_lim_w), abs(self.pack.cold.sum_of_q_lim_w)])  # abs, т.к. разные знаки
        q_act_w = q_lim_w * eps_hex

        if q_lim_w == self.pack.cold.sum_of_q_lim_w:  # self.pack.cold.sum_of_q_lim_w > 0
            pack_cold_has_min_q_w = True
            pack_active = self.pack.hot
            pack_passive = self.pack.cold
        else:
            pack_hot_has_min_q_w = True
            pack_active = self.pack.cold
            pack_passive = self.pack.hot

        if abs(self.pack.hot.sum_of_q_lim_w - self.pack.cold.sum_of_q_lim_w) < 0.0001:
            print(f'Heat exchanger hot and flow packs are perfectly balanced with heat capacities:')
            print(f'hot pack of  : {self.pack.hot.sum_of_q_lim_w}, [W]')
            print(f'cold pack of : {self.pack.cold.sum_of_q_lim_w}, [W]')
            sys.exit('case is not considered yet. program terminated in "_compare_packs__pack_flag_q_lim_w_q_act_w"')

        return q_lim_w, q_act_w, pack_active, pack_passive, pack_hot_has_min_q_w, pack_cold_has_min_q_w

    # def compare_packs_set_q_lim_w_q_act_w(self, eps_hex):
    #
    #     pack_active_pointer, pack_passive_pointer = self.set_active_passive_packs()

    def calc_qt(self, sorting_func) -> bool:
        # расчет кривых t(x), h(x), s(x) и dq(x) потоков внутри пассивного пакета
        # (у которого только один набор eps в eps_sets).
        self.pack.passive.thsdq_profiling(i_operating_mode=0)

        # программа сортировки sorting_func() - это внешняя функция, переданная в качестве аргумента. она возвращает
        # кортеж кортежей, где каждый внутренний кортеж, скажем для случая, когда параметр - производства энтропии - это
        # (i, ds), где i - индекс рабочего режима (индекс внутр. кортежа в self.pack.passive.flows_outlet_dat_in_sets),
        # для которого было рассчит. производство энтропии ds. Для случая ds внутренние кортежи отсортированы
        # с возрастанием ds. будем последовательно перебирать i из внутр. кортежей и проводить для них расчет qt
        operating_modes_sorted = sorting_func()

        qt_is_operable = False  # иначе return qt_is_operable говорит, что qt_is_operable не был присвоен
        for operating_mode in operating_modes_sorted:
            i_mode = operating_mode[0]

            # расчет кривых t(x), h(x), s(x) и dq(x) потоков внутри активного пакета
            # (у которого много наборов eps в eps_sets; eps_sets(i_mode)).
            self.pack.active.thsdq_profiling(i_operating_mode=i_mode)

            # calculate inter flows temperature differences as a matrix:
            qt_is_operable, \
                self.dt_inter_flows, \
                self.dt_min_inter_flows, \
                self.dt_max_inter_flows, \
                self.dt_avr_inter_flows = self.calc_dt_inter_flows(dt_min=self.dt_inter_flows_min)
            if qt_is_operable:
                # расчет производства энтропии по сечениям ТО. т.е. сумма произв. ds по всем потокам обоих пакетов
                self.ds_wk, self.ds_wk_sum = self.calc_ds_wk_distribution()
                print(f'for i_mode = {i_mode} qt-diagram is OPERABLE')
                break
            else:
                print(f'for i_mode = {i_mode} qt-diagram is INOPERABLE')

        return qt_is_operable

        #     if not no_t_crossing_in_pack:   # кривые t(x) потоков активного пакета пересеклись между собой: FAILURE
        #         qt_inf = 't_crossing_in_active_pack'
        #         qt_is_operable = False
        #     else:   # кривые t(x) потоков активного пакета не пересеклись между собой: OK
        #         # проверяем нет ли пересечения кривых t(x) самого холодного потока в "горячем" пакете
        #         # и самого горячего потока в "холодном" пакете: dt_min >= self.dt_inter_flows_min
        #         no_t_crossing_inter_pack, self.dt_inter_pack = self.dtmin_inter_packs_is_ok(dt_min=
        #                                                                                     self.dt_inter_flows_min)
        #         if not no_t_crossing_inter_pack:    # dt_min < self.dt_inter_flows_min
        #             qt_inf = 't_crossing_inter_pack'
        #             qt_is_operable = False
        #         else:   # для текущего i_mode рабочего режима qt-диаграмма полностью работоспособна
        #             break   # прекращаем перебор рабочих режимов
        #
        # # расчет производства энтропии по сечениям ТО. т.е. сумма произв. ds по всем потокам обоих пакетов
        # self.ds_wk, self.ds_wk_sum = self.calc_ds_wk_distribution()

        # print(f' qt_is_operable = {qt_is_operable}')
        # print(f' qt_inf = {qt_inf}')
        # sys.exit('ok')

        # print(f' ds_passive_pack = {pack_passive_pointer.ds_wk_sum}'
        #       f' ds_active_pack = {pack_active_pointer.ds_wk_sum}')
        # # f' ds_sum = {pack_active_pointer.ds_wk_sum+pack_passive_pointer.ds_wk_sum}')
        #
        # print(f'---------------------------------------------------------------------')
        # print(
        #     f'  i_oper_mode         eps_set       ds   t_out[0]  t_out[1]   dt_out[0]  dt_out[1]  ds[0]    ds[1]       qt_inf')
        # print(f'---------------------------------------------------------------------')
        # for i_oper_mode in range(len(pack_active_pointer.flows_outlet_dat_in_sets)):
        #
        #     qt_is_operable = True
        #     qt_inf = 'ok'
        #
        #     no_t_crossing_in_active_pack = pack_active_pointer.calc_thsdq_distribution_in_active_pack(i_calc=
        #                                                                                               i_oper_mode)
        #
        #     # no_t_crossing_in_active_pack: флаг наличия пересечения тем-ных крив. потоков в одном пакете
        #     if not no_t_crossing_in_active_pack:
        #         qt_inf = 't_crossing_in_active_pack'
        #         qt_is_operable = False
        #     # no_t_crossing_inter_pack: флаг наличия пересеч. тем-ных крив. самого холодного потока в "горячем" пакете
        #     # и самого горячего потока в "холодном" пакете
        #     else:
        #         no_t_crossing_inter_pack, self.dt_inter_pack = self.dtmin_inter_packs_is_ok(dt_min=
        #                                                                                     self.dt_inter_flows_min)
        #         if not no_t_crossing_inter_pack:
        #             qt_inf = 't_crossing_inter_pack'
        #             qt_is_operable = False

        # break

        # if qt_is_operable:
        #     break
        # печать рез-тов расч. температур потоков актив. пакета: i_sec t[K] dq[W] (для всех потоков) dq_sum (пакета)
        # pack_active_pointer.print_t_q_distribution()

        # for calculated qt-diagram calculate inter flows temperature differences as a matrix:
        # self.dt_inter_flows, self.dt_min_inter_flows, self.dt_max_inter_flows, self.dt_avr_inter_flows
        # self.dt_inter_flows, \
        # self.dt_min_inter_flows, \
        # self.dt_max_inter_flows, \
        # self.dt_avr_inter_flows = self.calc_dt_inter_flows()

        # ----------------------------------------------------------------------------------------------------------
        # этот фрагмент нужен для анализа результатов и является избыточным
        # eps_set = [round(flow.dq_w_sum / pack_active_pointer.dq_w_sum, 2) for flow in pack_active_pointer.flow]
        # dt_lim = [round(flow.outlet_lim.t_k - flow.t_k[-1], 2) for flow in pack_active_pointer.flow]
        # t_out = [round(flow.t_k[-1], 2) for flow in pack_active_pointer.flow]
        # ds = pack_active_pointer.flows_outlet_dat_in_sets[i_oper_mode][1]
        # ds_flow = [round((flow.s_jmolk[-1] - flow.s_jmolk[0]) * flow.m_mols, 2) for flow in
        #            pack_active_pointer.flow]
        # if qt_inf != 'ok':
        #     print(f'       {i_oper_mode:2d}'
        #           f'           {eps_set[0]:.2f}   {eps_set[1]:.2f}'
        #           f'{ds:8.2f} '
        #           f'  {t_out[0]:.2f}     {t_out[1]:.2f}'
        #           f'   {dt_lim[0]:7.2f}   {dt_lim[1]:7.2f}'
        #           f'   {ds_flow[0]:7.2f}   {ds_flow[1]:7.2f}'
        #           f'{qt_inf:>30s} ', end='')
        #     # print(f'
        #     # только для случая одного гор.потока и двух холодных
        #     i_hot = 0
        #     j_cold = 0
        #     print(f'\t {self.dt_min_inter_flows[i_hot][j_cold]:6.2f}'
        #           f'\t {self.dt_max_inter_flows[i_hot][j_cold]:6.2f}'
        #           f'\t {self.dt_avr_inter_flows[i_hot][j_cold]:6.2f}'
        #           f'\t {self.dt_min_inter_flows[i_hot][j_cold + 1]:6.2f}'
        #           f'\t {self.dt_max_inter_flows[i_hot][j_cold + 1]:6.2f}'
        #           f'\t {self.dt_avr_inter_flows[i_hot][j_cold + 1]:6.2f}')
        # ----------------------------------------------------------------------------------------------------------

        # i_hot = 0
        # j_cold = 0
        # if qt_is_operable:
        #     print(f' {i_calc:2d}'
        #           f'\t {pack_active_pointer.flows_outdat_in_sets[i_calc][1]:8.2f}'  # ds
        #           f'\t {self.dt_min_inter_flows[i_hot][j_cold]:6.2f}'
        #           f'\t {self.dt_max_inter_flows[i_hot][j_cold]:6.2f}'
        #           f'\t {self.dt_avr_inter_flows[i_hot][j_cold]:6.2f}'
        #           f'\t {self.dt_min_inter_flows[i_hot][j_cold + 1]:6.2f}'
        #           f'\t {self.dt_max_inter_flows[i_hot][j_cold + 1]:6.2f}'
        #           f'\t {self.dt_avr_inter_flows[i_hot][j_cold + 1]:6.2f}'
        #           f'\t {self.dt_min_inter_flows[i_hot + 1][j_cold]:6.2f}'
        #           f'\t {self.dt_max_inter_flows[i_hot + 1][j_cold]:6.2f}'
        #           f'\t {self.dt_avr_inter_flows[i_hot + 1][j_cold]:6.2f}'
        #           f'\t {self.dt_min_inter_flows[i_hot + 1][j_cold + 1]:6.2f}'
        #           f'\t {self.dt_max_inter_flows[i_hot + 1][j_cold + 1]:6.2f}'
        #           f'\t {self.dt_avr_inter_flows[i_hot + 1][j_cold + 1]:6.2f}'
        #           )
        # else:
        #     print(f' {i_calc:2d}'
        #           f'\t {pack_active_pointer.flows_outdat_in_sets[i_calc][1]:8.2f}'  # ds
        #           f'\t {self.dt_min_inter_flows[i_hot][j_cold]:6.2f}'
        #           f'\t {self.dt_max_inter_flows[i_hot][j_cold]:6.2f}'
        #           f'\t {self.dt_avr_inter_flows[i_hot][j_cold]:6.2f}'
        #           f'\t {self.dt_min_inter_flows[i_hot][j_cold + 1]:6.2f}'
        #           f'\t {self.dt_max_inter_flows[i_hot][j_cold + 1]:6.2f}'
        #           f'\t {self.dt_avr_inter_flows[i_hot][j_cold + 1]:6.2f}'
        #           f'\t {self.dt_min_inter_flows[i_hot + 1][j_cold]:6.2f}'
        #           f'\t {self.dt_max_inter_flows[i_hot + 1][j_cold]:6.2f}'
        #           f'\t {self.dt_avr_inter_flows[i_hot + 1][j_cold]:6.2f}'
        #           f'\t {self.dt_min_inter_flows[i_hot + 1][j_cold + 1]:6.2f}'
        #           f'\t {self.dt_max_inter_flows[i_hot + 1][j_cold + 1]:6.2f}'
        #           f'\t {self.dt_avr_inter_flows[i_hot + 1][j_cold + 1]:6.2f}'
        #           f'\t   qt_is_not_operable'
        #           )

    def calc_ds_wk_distribution(self) -> [[], float]:
        ds_wk = [0.0]
        for i in range(1, self.n_qt_sect + 1):
            ds_wk_section_hot_pack = 0
            for flow in self.pack.hot.flow:
                ds_wk_section_hot_pack += flow.ds_wk[i - 1]
            ds_wk_section_cold_pack = 0
            for flow in self.pack.cold.flow:
                ds_wk_section_cold_pack += flow.ds_wk[i]
            ds_wk.append(ds_wk_section_hot_pack + ds_wk_section_cold_pack)
        ds_wk_sum = sum(ds_wk)
        return ds_wk, ds_wk_sum

    def calc_dt_inter_flows(self, dt_min):
        # расчет разниц температур dt между горячими и холодными потоками
        # dt_inter_flows[i][j][k] - разница температур между i-ым гор.потоком и j-ым холодным в k-ом сечении потока
        # dt_min[i][j], _max[i][j] и _avr[i][j]: минимальная, максимальная и средняя разности тем-р между
        # i-ым гор.потоком и j-ым холодным
        qt_is_operable = True
        dt_inter_flows = np.empty((self.pack.hot.n_flows, self.pack.cold.n_flows, self.n_qt_nods))
        dt_min_inter_flows = np.empty((self.pack.hot.n_flows, self.pack.cold.n_flows))
        dt_max_inter_flows = np.empty((self.pack.hot.n_flows, self.pack.cold.n_flows))
        dt_avr_inter_flows = np.empty((self.pack.hot.n_flows, self.pack.cold.n_flows))

        for i_hot, flow_hot in enumerate(self.pack.hot.flow):
            for j_cold, flow_cold in enumerate(self.pack.cold.flow):
                dt_inter_flows[i_hot][j_cold] = flow_hot.t_k - flow_cold.t_k  # сразу по всем сечениям потоков
                dt_min_inter_flows[i_hot][j_cold] = min(dt_inter_flows[i_hot][j_cold])
                if dt_min_inter_flows[i_hot][j_cold] < dt_min:
                    qt_is_operable = False
                dt_max_inter_flows[i_hot][j_cold] = max(dt_inter_flows[i_hot][j_cold])
                dt_avr_inter_flows[i_hot][j_cold] = np.average(dt_inter_flows[i_hot][j_cold])

        # print(dt_inter_flows)
        # print(dt_min_inter_flows)
        # print(dt_max_inter_flows)
        # print(dt_avr_inter_flows)
        return qt_is_operable, dt_inter_flows, dt_min_inter_flows, dt_max_inter_flows, dt_avr_inter_flows

    def dtmin_inter_packs_is_ok(self, dt_min: float = 0.0) -> [bool, []]:
        # находим минимальную разность температур между самым холодным потоком в "горячем" пакете pack.hot.flow[-1] и
        # самым горячим в "холодном" пакете pack.cold.flow[0]. до этого момента уже надо было убедиться в том,
        # что потоки одного пакета не пересекаются между собой!!!
        # ВАЖНО: списки потоков для пакетов были отсортированы по убыванию t_in (от большей t_in к меньшей)
        # еще при их создании в ..._flows_inp
        is_ok = True
        dlt_t = self.pack.hot.flow[-1].t_k - self.pack.cold.flow[0].t_k  # self.flow[i].t_k is np.array

        # при eps_hex = 1 температуры потоков "пассивного" пакет на выходе должны быть равны предельной температуре
        # потоков "активного" пакета на входе.
        # тогда, при расчете min(dlt_t) не будем учитывать крайних точек:
        if self.eps_hex == 1.0:
            dt_min_inter_packs = min(dlt_t[1:-1])
        else:
            dt_min_inter_packs = min(dlt_t)

        if dt_min_inter_packs < dt_min:
            is_ok = False

        return is_ok, dlt_t

    def _set_flows_outdat_combinations_useful_for_operating_modes(self):
        return self.pack.active._combine_flows_outdat_in_sets(), self.pack.passive._combine_flows_outdat_in_sets()

    def split_into_sections_with_dt(self) -> [int, int]:
        _n_sections = int((self.pack.hot.t_in_lim_k - self.pack.cold.t_in_lim_k) / self.dt_sect)
        _n_elements = _n_sections + 1
        return _n_sections, _n_elements

    def _calc_flows_outlet_dat_with_eps(self):
        # каждый поток пакета имеет свой кортеж eps: (eps_0, eps_1,...,eps_N ). для каждого eps_i проводим расчет
        # соответствующих выходных данных и сохраняем их в виде кортежа данных: (dat_0, dat_1,...,dat_N )
        # кортежи данных "для потоков" сохраняем в общем для них кортеже "для пакета":
        #                       ( (dat_0, dat_1,...,dat_N), (dat_0, dat_1,...,dat_N ),...)
        # т.е. данные для i-го потока пакета и его j-го eps будем в дальнейшем находить по индексу dat[i][j]
        outlet_dat_pack_active = self.pack.active.calc_flows_out_thsds_with_eps_qact(q_act_w=self.q_act_w)
        outlet_dat_pack_passive = self.pack.passive.calc_flows_out_thsds_with_eps_qact(q_act_w=self.q_act_w)
        return outlet_dat_pack_active, outlet_dat_pack_passive

    def _set_sections_in_flows(self, n_sections: int, n_nodes: int) -> None:
        self.pack.hot._set_sections_in_flows(n_sections=n_sections, n_nodes=n_nodes)
        self.pack.cold._set_sections_in_flows(n_sections=n_sections, n_nodes=n_nodes)
