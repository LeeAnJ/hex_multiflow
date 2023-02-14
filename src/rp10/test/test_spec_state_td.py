from prettytable import PrettyTable
import data.butane_ethane_methane_spec_state as sample
from calc_td import calc_td_and_compare_with_sample

def test_spec_state_td(mixture):

    # print output data into prettytable: set table instance and table header
    table_1 = PrettyTable(padding_width=2)
    col1_header = "Spec.-state f(t,d,x): butane-ethane-methane"
    col2_header = "EPSmin, %"
    col3_header = "EPSmax, %"
    table_1.field_names = [col1_header, col2_header, col3_header]
    table_1_rows = list()

    # SUBCOOLED LIQUID PHASE -------------------------------------------------------------------------------------------
    # get mixture temperature from sample dict: sample.rp10_spec_state_liq
    t_units = 'k'
    t = sample.rp10_spec_state_liq['t'][t_units]
    # get mixture density from sample dict: sample.rp10_spec_state_liq
    d_units = 'moll'
    d = sample.rp10_spec_state_liq['d'][d_units]
    x_units = 'molmol'
    x = sample.rp10_spec_state_liq['x'][x_units]

    # calc. f(t,d,x) for subcooled liquid phase and compare results with sample
    data_list = calc_td_and_compare_with_sample(fluid=mixture,
                                                t=t, t_units=t_units,
                                                d=d, d_units=d_units,
                                                x=x, x_units=x_units)
    # results: compare blk with sample
    for list_ in data_list: table_1_rows.append(list_)

    # 2PH --------------------------------------------------------------------------------------------------------------
    # input data from blk_sample
    t = sample.rp10_spec_state_2h['t'][t_units]
    d = sample.rp10_spec_state_2h['d'][d_units]
    data_list = calc_td_and_compare_with_sample(fluid=mixture,
                                                t=t, t_units=t_units,
                                                d=d, d_units=d_units,
                                                x=x, x_units=x_units)
    # results: compare blk, blk_liq and blk_vap with samples
    for list_ in data_list: table_1_rows.append(list_)

    # SUPER HEATED VAPOR -----------------------------------------------------------------------------------------------
    # input data from blk_sample
    t = sample.rp10_spec_state_vap['t'][t_units]
    d = sample.rp10_spec_state_vap['d'][d_units]
    data_list = calc_td_and_compare_with_sample(fluid=mixture,
                                                t=t, t_units=t_units,
                                                d=d, d_units=d_units,
                                                x=x, x_units=x_units)
    # results: compare blk with sample
    for list_ in data_list: table_1_rows.append(list_)

    # 2PH --------------------------------------------------------------------------------------------------------------
    # input data from blk_sample in USER units
    t_units = 'c'
    d_units = 'kgm3'
    x_units = 'kgkg'
    t = sample.rp10_spec_state_2h['t'][t_units]
    d = sample.rp10_spec_state_2h['d'][d_units]
    x = sample.rp10_spec_state_2h['x'][x_units]
    data_list = calc_td_and_compare_with_sample(fluid=mixture,
                                                t=t, t_units=t_units,
                                                d=d, d_units=d_units,
                                                x=x, x_units=x_units)
    # results: compare blk, blk_liq and blk_vap with samples
    for list_ in data_list: table_1_rows.append(list_)

    table_1.add_rows(table_1_rows)

    table_1.align[col1_header] = "l"
    table_1.align[col2_header] = "l"
    table_1.align[col3_header] = "l"

    print(table_1)
