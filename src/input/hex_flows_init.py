import sys

from src.hex.single_flow.flow_abc import set_isobar_for_hex_flow, FlowInputData
from src.hex.single_flow.flow_hot import HexFlowHot
from src.hex.single_flow.flow_cold import HexFlowCold
from src.input.hot_flows_inp import flow_in_hot_pack
from src.input.cold_flows_inp import flow_in_cold_pack
from copy import copy

# Пакет горячих потоков:
n_hot_flows = len(flow_in_hot_pack)
#   предполагается, что в ТО горячие потоки могут охладиться минимум до t_flows_cold_in_max,
#   а холодные потоки нагреться максимум до t_flows_hot_in_min (тогда не будет пересечения кривых t(x)

# устаревший код для расчета минимальной температуры горячих потоков на входе
# t_flows_hot_in_min = min([flow_in_hot_pack[i]['t_inlet'].value for i in range(n_hot_flows)])

# т.к. потоки в "горячем" пакете уже отсортированы по убыванию тем-ры на входе: см. hot_flows_inp.py
# то мин. тем-ра - это температура последнего потока в списке flow_in_hot_pack
t_flows_hot_in_min = copy(flow_in_hot_pack[n_hot_flows - 1]['t_inlet'])

# Пакет холодных потоков:
n_cold_flows = len(flow_in_cold_pack)

# устаревший код для расчета максимальная температура холодных потоков на входе
# t_flows_cold_in_max = max([flow_in_cold_pack[i]['t_inlet'].value for i in range(n_cold_flows)])

# т.к. потоки в "холодном" пакете уже отсортированы по убыванию тем-ры на входе: см. cold_flows_inp.py
# то максимальная тем-ра - это температура первого потока в списке flow_in_cold_pack
t_flows_cold_in_max = copy(flow_in_cold_pack[0]['t_inlet'])

# список вход. параметров для всех гор. потоков c dataclass FlowInputData :
#     flow_hot = [ ((isobar,m,t,p), t_out_lim), ((isobar,m,t,p), t_out_lim), ... ]
flow_hot = []
for i in range(n_hot_flows):
    #   для построения изобары использую ф-цию set_isobar_for_hex_flow из flow_abc.py, а не просто Isobar,
    #   т.к. крайние тем-ры t_min, t_max надо "сдвинуть" на dt_end, а то в них не будет работать лин. интерполяция
    #   в конструктор изобары Isobar можно передавать p и t как они были заданы в hot/cold_flows_inp, т.е. тут
    #   допустима любая размерность p и t, т.к. она будет приведена к internal units в конструкторе изобары Isobar
    flow_hot.append(HexFlowHot(FlowInputData(isobar=set_isobar_for_hex_flow(fluid=flow_in_hot_pack[i]['fluid'],
                                                                            p=flow_in_hot_pack[i]['p'],
                                                                            t_min=t_flows_cold_in_max,
                                                                            t_max=flow_in_hot_pack[i]['t_inlet']),
                                             mass_rate=copy(flow_in_hot_pack[i]['mass_rate']),
                                             t_in=copy(flow_in_hot_pack[i]['t_inlet']),
                                             p_in=copy(flow_in_hot_pack[i]['p'])),
                               t_out_lim=t_flows_cold_in_max))

# список вход. параметров для всех холодных потоков c dataclass FlowInputData :
#     flow_cold = [ ((isobar,m,t,p), t_out_lim), ((isobar,m,t,p), t_out_lim), ... ]
flow_cold = []
for i in range(n_cold_flows):
    #   для построения изобары использую ф-цию set_isobar_for_hex_flow из flow_abc.py, а не просто Isobar,
    #   т.к. крайние тем-ры t_min, t_max надо "сдвинуть" на dt_end, а то в них не будет работать лин. интерполяция
    #   в конструктор изобары Isobar можно передавать p и t как они были заданы в hot/cold_flows_inp, т.е. тут
    #   допустима любая размерность p и t, т.к. она будет приведена к internal units в конструкторе изобары Isobar
    flow_cold.append(HexFlowCold(FlowInputData(isobar=set_isobar_for_hex_flow(fluid=flow_in_cold_pack[i]['fluid'],
                                                                              p=flow_in_cold_pack[i]['p'],
                                                                              t_min=flow_in_cold_pack[i]['t_inlet'],
                                                                              t_max=t_flows_hot_in_min),
                                               mass_rate=copy(flow_in_cold_pack[i]['mass_rate']),
                                               t_in=copy(flow_in_cold_pack[i]['t_inlet']),
                                               p_in=copy(flow_in_cold_pack[i]['p'])),
                                 t_out_lim=t_flows_hot_in_min))

# sys.exit('ok')
