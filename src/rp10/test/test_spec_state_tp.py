from prettytable import PrettyTable
import data.butane_ethane_methane_spec_state as sample
from calc_tp import calc_tp_and_compare_with_sample

def test_spec_state_tp(mixture):

    # print output data into prettytable: set table instance and table header
    table_1 = PrettyTable(padding_width=2)
    col1_header = "Spec.-state f(t,p,x): butane-ethane-methane"
    col2_header = "EPSmin, %"
    col3_header = "EPSmax, %"
    table_1.field_names = [col1_header, col2_header, col3_header]

    table_1_rows = list()

    # get mixture temperature from sample dict: sample.rp10_spec_state_liq
    t_units = 'k'
    t = sample.rp10_spec_state_liq['t'][t_units]
    # get mixture pressure from sample dict: sample.rp10_spec_state_liq
    p_units = 'kpa'
    p = sample.rp10_spec_state_liq['p'][p_units]
    x_units = 'molmol'
    x = sample.rp10_spec_state_liq['x'][x_units]

    # calc. f(t,p,x) for subcooled liquid phase and compare results with sample
    data_list = calc_tp_and_compare_with_sample(fluid=mixture, t=t, t_units=t_units, p=p, p_units=p_units, x=x, x_units=x_units)
    for list_ in data_list: table_1_rows.append(list_)

    t = sample.rp10_spec_state_2h['t'][t_units]
    data_list = calc_tp_and_compare_with_sample(fluid=mixture, t=t, t_units=t_units, p=p, p_units=p_units, x=x, x_units=x_units)
    for list_ in data_list: table_1_rows.append(list_)

    t = sample.rp10_spec_state_vap['t'][t_units]
    data_list = calc_tp_and_compare_with_sample(fluid=mixture, t=t, t_units=t_units, p=p, p_units=p_units, x=x, x_units=x_units)
    for list_ in data_list: table_1_rows.append(list_)

    t_units = 'c'
    t = sample.rp10_spec_state_liq['t'][t_units]
    # get mixture pressure from sample dict: sample.rp10_spec_state_liq
    p_units = 'bar'
    p = sample.rp10_spec_state_liq['p'][p_units]
    x_units = 'kgkg'
    x = sample.rp10_spec_state_liq['x'][x_units]
    data_list = calc_tp_and_compare_with_sample(fluid=mixture, t=t, t_units=t_units, p=p, p_units=p_units, x=x, x_units=x_units)
    for list_ in data_list: table_1_rows.append(list_)

    table_1.add_rows(table_1_rows)

    table_1.align[col1_header] = "l"
    table_1.align[col2_header] = "l"
    table_1.align[col3_header] = "l"

    print(table_1)

    # for row in table_1_rows:
    #     print(row)
