import sys

from src.input.hex_flows_init import flow_hot, flow_cold
from src.output.table.p_table import print_input_data_as_table, print_t_dq_with_cross_section_as_table, \
    print_dt_with_cross_section_as_table, print_t_ds_with_cross_section_as_table
from src.output.plot_fig.plot_fig import plot_t_with_cross_sections
from src.hex.hex import Hex, HexFlows, HexData

# список экз.!!! горячих потоков flow_hot "собран" и отсортирован по убыванию t_in в hot_flows_inp
# список экз.!!! холодных потоков flow_cold "собран" и отсортирован по убыванию t_in в cold_flows_inp
rhex = Hex(flows_sorted_with_tin=HexFlows(hot=flow_hot, cold=flow_cold),
           data=HexData(dt_as=10.0,  # шаг по тем-ре вдоль потока (=> кол-во секций ТО)
                        dt_cs=3.0,  # минимально допустимая разность тем-р между хол. и гор. потоками
                        eps_hex=0.950,  # эффектив. ТО (потери): eps_hex = Q_act/Q_max, Q_max=min(Q_sum_hot,Q_sum_cold)
                        dlt_eps=0.01))  # шаг (точность) по eps для потока: eps[i][j]=Q[i]/Q_act

# qt-диаграмма ТО: расчет
# если сорт.ф-ция rhex.calc_ds_with_operating_modes, то рабочие режимы отсортированы по возраст. энтропии,сгенерир. в ТО
qt_is_operable = rhex.calc_qt(sorting_func=rhex.calc_ds_with_operating_modes)

# печать результатов в графическом виде
plot_t_with_cross_sections(hex=rhex)
# sys.exit('ok')

# печать результатов в табличном виде
print_input_data_as_table(hex=rhex)  # Table: Hex input data: fluids, mass rates,
print_t_dq_with_cross_section_as_table(hex=rhex,
                                       qt_is_operable_flag=qt_is_operable)  # Table: t,dq = f(index_cross_section)
print_dt_with_cross_section_as_table(hex=rhex)  # Table: dt = f(index_cross_section)
# print_t_ds_with_cross_section_as_table(hex=rhex)                          # Table: t,ds = f(index_cross_section)
