from test_functions import compare_data_dicts
import data.butane_ethane_methane_spec_state as sample


def compare_calculated_and_sample_data(fluid=None,
                                       x_name=None, x_units=None,
                                       y_name=None, y_units=None,
                                       z_name=None, z_units=None):
    func_name = " spec_state("+x_name+"["+x_units+"], "+y_name+"["+y_units+"], "+z_name+"["+z_units+"])"
    list_of_results = list()

    if fluid.error.index > 0:
        fluid.error.print_and_terminate()
    else:
        fluid.convert_dataset_to_user_units(flag_data='blk')
        phase_symbol = fluid.state.get_phase_symbol()
        if phase_symbol == 'l':
            list_ = compare_data_dicts(dict_model=sample.rp10_spec_state_liq, dict_calc=fluid.state.data['blk'])
            list_.insert(0, 'l-phase  '+func_name+'['+'blk'+']       :')
            list_of_results.append(list_)
        elif phase_symbol == 'v':
            list_ = compare_data_dicts(dict_model=sample.rp10_spec_state_vap, dict_calc=fluid.state.data['blk'])
            list_.insert(0, 'v-phase  '+func_name+'['+'blk'+']       :')
            list_of_results.append(list_)
        else:   # phase_symbol = 'l_v'
            # blk
            list_ = compare_data_dicts(dict_model=sample.rp10_spec_state_2h, dict_calc=fluid.state.data['blk'])
            list_.insert(0, '2ph.     '+func_name+'[' + 'blk' + ']       :')
            list_of_results.append(list_)
            # blk_liq
            fluid.convert_dataset_to_user_units(flag_data='blk_liq')
            list_ = compare_data_dicts(dict_model=sample.rp10_spec_state_2h_liq, dict_calc=fluid.state.data['blk_liq'])
            list_.insert(0, '2ph_liq. '+func_name+'[' + 'blk_liq' + ']   :')
            list_of_results.append(list_)
            # blk_vap
            fluid.convert_dataset_to_user_units(flag_data='blk_vap')
            list_ = compare_data_dicts(dict_model=sample.rp10_spec_state_2h_vap, dict_calc=fluid.state.data['blk_vap'])
            list_.insert(0, '2ph_vap. '+func_name+'[' + 'blk_vap' + ']   :')
            list_of_results.append(list_)

    return list_of_results
