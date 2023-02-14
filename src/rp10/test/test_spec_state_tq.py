from prettytable import PrettyTable
import data.butane_ethane_methane_spec_state as sample
from calc_tq import calc_tq_and_compare_with_sample


def test_spec_state_tq(mixture):

    # print output data into prettytable: set table instance and table header
    table_1 = PrettyTable(padding_width=2)
    col1_header = "Spec.-state f(t,q,x): butane-ethane-methane"
    col2_header = "EPSmin, %"
    col3_header = "EPSmax, %"
    table_1.field_names = [col1_header, col2_header, col3_header]
    table_1_rows = list()


    # 2PH --------------------------------------------------------------------------------------------------------------
    # input data from blk_sample
    t_units = 'k'
    t = sample.rp10_spec_state_2h['t'][t_units]
    # get mixture quality from sample dict: sample.rp10_spec_state_2ph
    q_units = 'molmol'
    q = sample.rp10_spec_state_2h['q'][q_units]
    x_units = 'molmol'
    x = sample.rp10_spec_state_2h['x'][x_units]
    data_list = calc_tq_and_compare_with_sample(fluid=mixture,
                                                t=t, t_units=t_units,
                                                q=q, q_units=q_units,
                                                x=x, x_units=x_units)
    # results: compare blk, blk_liq and blk_vap with samples
    for list_ in data_list:
        table_1_rows.append(list_)

    # 2PH --------------------------------------------------------------------------------------------------------------
    # input data from blk_sample in USER units
    t_units = 'c'
    q_units = 'kgkg'
    x_units = 'kgkg'
    t = sample.rp10_spec_state_2h['t'][t_units]
    q = sample.rp10_spec_state_2h['q'][q_units]
    x = sample.rp10_spec_state_2h['x'][x_units]
    data_list = calc_tq_and_compare_with_sample(fluid=mixture,
                                                t=t, t_units=t_units,
                                                q=q, q_units=q_units,
                                                x=x, x_units=x_units)
    # results: compare blk, blk_liq and blk_vap with samples
    for list_ in data_list:
        table_1_rows.append(list_)

    table_1.add_rows(table_1_rows)

    table_1.align[col1_header] = "l"
    table_1.align[col2_header] = "l"
    table_1.align[col3_header] = "l"

    print(table_1)
