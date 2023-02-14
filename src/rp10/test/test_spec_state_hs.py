from prettytable import PrettyTable
import data.butane_ethane_methane_spec_state as sample
from calc_hs import calc_hs_and_compare_with_sample


def test_spec_state_hs(mixture):

    # print output data into prettytable: set table instance and table header
    table_1 = PrettyTable(padding_width=2)
    col1_header = "Spec.-state f(h,s,x): butane-ethane-methane"
    col2_header = "EPSmin, %"
    col3_header = "EPSmax, %"
    table_1.field_names = [col1_header, col2_header, col3_header]
    table_1_rows = list()

    # SUBCOOLED LIQUID PHASE -------------------------------------------------------------------------------------------
    # get mixture temperature from sample dict: sample.rp10_spec_state_liq
    h_units = 'jmol'
    h = sample.rp10_spec_state_liq['h'][h_units]
    # get mixture entropy from sample dict: sample.rp10_spec_state_liq
    s_units = 'jmolk'
    s = sample.rp10_spec_state_liq['s'][s_units]
    x_units = 'molmol'
    x = sample.rp10_spec_state_liq['x'][x_units]

    # calc. f(d,s,x) for subcooled liquid phase and compare results with sample
    data_list = calc_hs_and_compare_with_sample(fluid=mixture,
                                                h=h, h_units=h_units,
                                                s=s, s_units=s_units,
                                                x=x, x_units=x_units)
    # results: compare blk with sample
    for list_ in data_list:
        table_1_rows.append(list_)

    # 2PH --------------------------------------------------------------------------------------------------------------
    # input data from blk_sample
    # h = sample.rp10_spec_state_2h['h'][h_units]
    # s = sample.rp10_spec_state_2h['s'][s_units]
    # data_list = calc_hs_and_compare_with_sample(fluid=mixture,
    #                                             h=h, h_units=h_units,
    #                                             s=s, s_units=s_units,
    #                                             x=x, x_units=x_units)
    # # results: compare blk, blk_liq and blk_vap with samples
    # for list_ in data_list:
    #     table_1_rows.append(list_)

    # SUPER HEATED VAPOR -----------------------------------------------------------------------------------------------
    # input data from blk_sample
    h = sample.rp10_spec_state_vap['h'][h_units]
    s = sample.rp10_spec_state_vap['s'][s_units]
    data_list = calc_hs_and_compare_with_sample(fluid=mixture,
                                                h=h, h_units=h_units,
                                                s=s, s_units=s_units,
                                                x=x, x_units=x_units)
    # results: compare blk with sample
    for list_ in data_list:
        table_1_rows.append(list_)

    # 2PH --------------------------------------------------------------------------------------------------------------
    # input data from blk_sample in USER units
    # h_units = 'jkg'
    # s_units = 'jkgk'
    # x_units = 'kgkg'
    # h = sample.rp10_spec_state_2h['h'][h_units]
    # s = sample.rp10_spec_state_2h['s'][s_units]
    # x = sample.rp10_spec_state_2h['x'][x_units]
    # data_list = calc_hs_and_compare_with_sample(fluid=mixture,
    #                                             h=h, h_units=h_units,
    #                                             s=s, s_units=s_units,
    #                                             x=x, x_units=x_units)
    # # results: compare blk, blk_liq and blk_vap with samples
    # for list_ in data_list:
    #     table_1_rows.append(list_)

    table_1.add_rows(table_1_rows)

    table_1.align[col1_header] = "l"
    table_1.align[col2_header] = "l"
    table_1.align[col3_header] = "l"

    print(table_1)
