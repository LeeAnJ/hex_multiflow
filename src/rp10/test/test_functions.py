import sys
from src.rp10.fluid.fluid_properties_class import RP10FluidData


def compare_data_dicts(dict_model: dict = None, dict_calc: dict = None) -> list:
    list_of_symbols = list()
    list_of_values = list()
    for _x_symbol, _x_dict in dict_model.items():
        for _x_units, _x_value in _x_dict.items():
            # if _x_symbol in ['x', 'q', 'mm']:
            if _x_symbol in ['x', 'q']:  # x- array, q=0 at sat.liq.curve
                pass
            else:
                if dict_model[_x_symbol][_x_units] is not None:
                    relative_error = (
                        abs(abs(dict_model[_x_symbol][_x_units] - dict_calc[_x_symbol][_x_units]) /
                            dict_model[_x_symbol][_x_units] * 100.0)
                    )
                    row_str = _x_symbol + ', [' + _x_units + ']'
                    list_of_symbols.append(row_str)
                    list_of_values.append(relative_error)
                    # print(row_str, dict_calc[_x_symbol][_x_units], dict_model[_x_symbol][_x_units], relative_error)
    # print(' ')
    relative_error_max = max(list_of_values)
    relative_error_max_index = list_of_values.index(relative_error_max)
    table_element_max = list_of_symbols[relative_error_max_index] + ':  ' + '{:.6f}'.format(relative_error_max)

    relative_error_min = min(list_of_values)
    relative_error_min_index = list_of_values.index(relative_error_min)
    table_element_min = list_of_symbols[relative_error_min_index] + ':  ' + '{:.6f}'.format(relative_error_min)

    # print(list_of_symbols)
    # print(list_of_values)

    return [table_element_min, table_element_max]


def print_sat_state_data(flag, data):  # flag in ['dew', 'bubble']
    # props_symbol = ['t', 'p', 'd', 'h', 'e', 's', 'cp', 'cv', 'q', 'w', 'eta', 'tcx', 'mm', 'x']
    # units_internal = ['k', 'kpa', 'moll', 'jmol', 'jmol', 'jmolk', 'jmolk', 'jmolk', 'molmol', 'ms', 'upas', 'wmk',
    #                   'gmol', 'molmol']
    # units_user = ['c', 'bar', 'kgm3', 'jkg', 'jkg', 'jkgk', 'jkgk', 'jkgk', 'kgkg', 'ms', 'upas', 'wmk',
    #               'gmol', 'kgkg']
    if flag not in ['dew', 'bubble']:
        sys.exit('input arg "flag" is out of rang in "set_sat_data_to_none" function of RP10FluidData class')

    for i in range(len(RP10FluidData.props_symbol)):
        symbol = RP10FluidData.props_symbol[i]
        units_internal = RP10FluidData.units_internal[i]
        units_user = RP10FluidData.units_user[i]
        print(symbol, data[flag][symbol][units_internal], ', [', units_internal, ']',
              data[flag][symbol][units_user], ', [', units_user, ']')
