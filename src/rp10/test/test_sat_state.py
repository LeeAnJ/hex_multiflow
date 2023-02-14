import sys

from src.rp10.fluid.fluid_class import RP10Fluid
import data.butane_ethane_methane_sat_state as mixture_sat_model
import data.butane_sat_state as butane_model
from prettytable import PrettyTable
from test_functions import compare_data_dicts


def butane_with_t(butane: RP10Fluid, sat_line: str, t_units: str) -> list:  # sat_line = 'v' or 'l'

    if sat_line == 'l':
        sat_str = 'liq'
    elif sat_line == 'v':
        sat_str = 'vap'
    else:
        sys.exit('wrong sat_line argument in func. test_butane_satt_state')

    if t_units == 'k':
        t_k__butane = butane_model.satt_liq_model_pattern['t']['k']
        butane.calc_sat_state(sat_curve_flag=sat_line, t=(t_k__butane, 'K'))
        line_tug = 'sat_state('+sat_str+', t [K])'
    elif t_units == 'c':
        t_c__butane = butane_model.satt_liq_model_pattern['t']['c']
        butane.calc_sat_state(sat_curve_flag=sat_line, t=(t_c__butane, 'c'))
        line_tug = 'sat_state('+sat_str+', t [oC])'
    else:
        sys.exit('wrong t_units argument in func. test_butane_satt_state')

    if butane.error.index > 0:
        butane.error.print_and_terminate()
    else:
        butane.convert_dataset_to_user_units(flag_data=sat_line)

        if sat_line == 'l':
            list_1 = compare_data_dicts(dict_model=butane_model.satt_liq_model_pattern,
                                        dict_calc=butane.state.data['bubble'])
        else:
            list_1 = compare_data_dicts(dict_model=butane_model.satt_vap_model_pattern,
                                        dict_calc=butane.state.data['dew'])
        list_1.insert(0, line_tug)
        return list_1


def butane_with_p(butane: RP10Fluid, sat_line: str, p_units: str) -> list:  # sat_line = 'v' or 'l'

    if sat_line == 'l':
        sat_str = 'liq'
    elif sat_line == 'v':
        sat_str = 'vap'
    else:
        sys.exit('wrong sat_line argument in func. test_butane_satp_state')

    if p_units == 'kpa':
        p_kpa_butane = butane_model.satp_liq_model_pattern['p']['kpa']
        butane.calc_sat_state(sat_curve_flag=sat_line, p=(p_kpa_butane, 'kPa'))
        line_tug = 'sat_state('+sat_str+', p [kPa])'
    elif p_units == 'bar':
        p_bar_butane = butane_model.satp_liq_model_pattern['p']['bar']
        butane.calc_sat_state(sat_curve_flag=sat_line, p=(p_bar_butane, 'bar'))
        line_tug = 'sat_state('+sat_str+', p [bar])'
    else:
        sys.exit('wrong p_units argument in func. test_butane_satp_state')

    if butane.error.index > 0:
        butane.error.print_and_terminate()
    else:
        butane.convert_dataset_to_user_units(flag_data=sat_line)

        if sat_line == 'l':
            list_1 = compare_data_dicts(dict_model=butane_model.satp_liq_model_pattern,
                                        dict_calc=butane.state.data['bubble'])
        else:
            list_1 = compare_data_dicts(dict_model=butane_model.satp_vap_model_pattern,
                                        dict_calc=butane.state.data['dew'])
        list_1.insert(0, line_tug)
        return list_1


def mixture_with_t(fluid: RP10Fluid, sat_line: str, t_units: str) -> list:  # sat_line = 'v' or 'l'
    if sat_line == 'l':
        sat_str = 'liq'
    elif sat_line == 'v':
        sat_str = 'vap'
    else:
        sys.exit('wrong sat_line argument in func. test_mixture_satt_state')

    if t_units == 'k':  # internal units: t, [K]; x, [mol/mol]
        t_k = mixture_sat_model.satt_liq_model_pattern['t']['k']
        fluid.calc_sat_state(sat_curve_flag=sat_line, t=(t_k, 'K'), x=(fluid.composition_molmol, 'molmol'))
        line_tug = sat_str+' : f( t[K],   x[mol/mol] )'
    elif t_units == 'c':  # user units: t, [oC]; x, [kg/kg]
        t_c = mixture_sat_model.satt_liq_model_pattern['t']['c']
        fluid.calc_sat_state(sat_curve_flag=sat_line, t=(t_c, 'C'), x=(fluid.composition_kgkg, 'kgkg'))
        line_tug = sat_str+' : f( t[C],   x[kg/kg] )'
    else:
        sys.exit('wrong t_units argument in func. test_mixture_satt_state')

    if fluid.error.index > 0:
        fluid.error.print_and_terminate()
    else:
        fluid.convert_dataset_to_user_units(flag_data=sat_line)
        if sat_line == 'l':
            list_1 = compare_data_dicts(dict_model=mixture_sat_model.satt_liq_model_pattern,
                                        dict_calc=fluid.state.data['bubble'])
        else:
            list_1 = compare_data_dicts(dict_model=mixture_sat_model.satt_vap_model_pattern,
                                        dict_calc=fluid.state.data['dew'])
        list_1.insert(0, line_tug)
        return list_1


def mixture_with_p(fluid: RP10Fluid, sat_line: str, p_units: str) -> list:  # sat_line = 'v' or 'l'
    # alternative choice: liquid (bubble) or vapour (dew) saturation line
    if sat_line == 'l':
        sat_str = 'liq'
    elif sat_line == 'v':
        sat_str = 'vap'
    else:
        sys.exit('wrong sat_line argument in func. test_mixture_satp_state')

    # alternative choice: calc. with internal or user units
    if p_units == 'kpa':  # internal units: p, [kPa]; x, [mol/mol]
        p_kpa = mixture_sat_model.satp_liq_model_pattern['p']['kpa']
        fluid.calc_sat_state(sat_curve_flag=sat_line, p=(p_kpa, 'kpa'), x=(fluid.composition_molmol, 'molmol'))
        line_tug = sat_str+' : f( p[kPa], x[mol/mol] )'
    elif p_units == 'bar':  # user units: p, [bar]; x, [kg/kg]
        p_bar = mixture_sat_model.satp_liq_model_pattern['p']['bar']
        fluid.calc_sat_state(sat_curve_flag=sat_line, p=(p_bar, 'bar'), x=(fluid.composition_kgkg, 'kgkg'))
        line_tug = sat_str+' : f( p[bar], x[kg/kg] )'
    else:
        sys.exit('wrong p_units argument in func. test_mixture_satp_state')

    # calculated and model sat. properties comparison
    if fluid.error.index > 0:
        fluid.error.print_and_terminate()
    else:
        fluid.convert_dataset_to_user_units(flag_data=sat_line)
        if sat_line == 'l':
            list_1 = compare_data_dicts(dict_model=mixture_sat_model.satp_liq_model_pattern,
                                        dict_calc=fluid.state.data['bubble'])
        else:
            list_1 = compare_data_dicts(dict_model=mixture_sat_model.satp_vap_model_pattern,
                                        dict_calc=fluid.state.data['dew'])
        list_1.insert(0, line_tug)
        return list_1
