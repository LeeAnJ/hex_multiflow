# import sys

from src.rp10.fluid.fluid_class import RP10Fluid
import src.rp10.test.test_sat_state as test_sat_state
import src.rp10.test.data.butane_ethane_methane_sat_state as mixture_model
from prettytable import PrettyTable

from src.rp10.test.test_spec_state_hs import test_spec_state_hs
from src.rp10.test.test_spec_state_pd import test_spec_state_pd
from src.rp10.test.test_spec_state_ph import test_spec_state_ph
from src.rp10.test.test_spec_state_ps import test_spec_state_ps
from src.rp10.test.test_spec_state_td import test_spec_state_td
from src.rp10.test.test_spec_state_tp import test_spec_state_tp
from src.rp10.test.test_spec_state_ds import test_spec_state_ds
from src.rp10.test.test_spec_state_tq import test_spec_state_tq

# --------------------------------------------set-BUTANE----------------------------------------------------------------
# initiate RP10Fluid pure fluid to be tested: ("butane",)
butane = RP10Fluid(names=("butane",))
# Output for butane as a prettytable: set table instance and table header.
table_butane = PrettyTable(padding_width=2)
table_butane.field_names = ["Sat. state at liq/vap curve",
                            "prop, [units] with : EPSmin, %",
                            "prop, [units] with : EPSmax, %"]

table_butane_rows = list()

# ---------------------------set-ZEOTROPIC MIXTURE of butane, ethane and methane----------------------------------------
x_units = 'molmol'
# mixture's component composition is taken from "model" data file: data.butane_ethane_methane_sat_state
# ("butane", "ethane", "methane"), (0.6, 0.3, 0.1), [mol/mol]
x = mixture_model.satt_liq_model_pattern['x'][x_units]
mixture = RP10Fluid(names=("butane", "ethane", "methane"), composition=(x, x_units))
# Output for mixture as a prettytable: set table instance and table header.
table_mixture = PrettyTable(padding_width=2)
table_mixture.field_names = ["Sat. state at liq/vap curve",
                             "prop, [units] with : EPSmin, %",
                             "prop, [units] with : EPSmax, %"]
table_mixture_rows = list()

# START SAT_STATES TESTING ---------------------------------------------------------------------------------------------
# calculate sat. properties for Butane and compare with model:
#   f(t): internal units t, [K]; p, [kPa]; user units t, [oC]; p, [bar];
table_butane_rows.append(test_sat_state.butane_with_t(butane=butane, sat_line='l', t_units='k'))
table_butane_rows.append(test_sat_state.butane_with_t(butane=butane, sat_line='l', t_units='c'))
table_butane_rows.append(test_sat_state.butane_with_t(butane=butane, sat_line='v', t_units='k'))
table_butane_rows.append(test_sat_state.butane_with_t(butane=butane, sat_line='v', t_units='c'))
#   f(p): internal units t, [K]; p, [kPa]; user units t, [oC]; p, [bar];
table_butane_rows.append(test_sat_state.butane_with_p(butane=butane, sat_line='l', p_units='kpa'))
table_butane_rows.append(test_sat_state.butane_with_p(butane=butane, sat_line='l', p_units='bar'))
table_butane_rows.append(test_sat_state.butane_with_p(butane=butane, sat_line='v', p_units='kpa'))
table_butane_rows.append(test_sat_state.butane_with_p(butane=butane, sat_line='v', p_units='bar'))

# Butane: output eps_min, eps_max as a pretty table
table_butane.add_rows(table_butane_rows)

table_butane.align["Sat. state at liq/vap curve"] = "l"
table_butane.align["prop, [units] with : EPSmin, %"] = "c"
table_butane.align["prop, [units] with : EPSmax, %"] = "c"

print(' ')
print(table_butane.get_string(title='Saturation state: butane  '
                                    '/ EPS = [ |prop.model[i] - prop.calc[i]| / prop.model[i]*100 ], % /'))

# calculate sat. properties for Mixture and compare with model:---------------------------------------------------------
#   f(t): internal units t, [K]; p, [kPa]; x, [mol/mol]; user units t, [oC]; p, [bar]; x, [kg/kg]
table_mixture_rows.append(test_sat_state.mixture_with_t(fluid=mixture, sat_line='l', t_units='k'))
table_mixture_rows.append(test_sat_state.mixture_with_t(fluid=mixture, sat_line='l', t_units='c'))
table_mixture_rows.append(test_sat_state.mixture_with_t(fluid=mixture, sat_line='v', t_units='k'))
table_mixture_rows.append(test_sat_state.mixture_with_t(fluid=mixture, sat_line='v', t_units='c'))
#   f(p): internal units t, [K]; p, [kPa]; x, [mol/mol]; user units t, [oC]; p, [bar]; x, [kg/kg]
table_mixture_rows.append(test_sat_state.mixture_with_p(fluid=mixture, sat_line='l', p_units='kpa'))
table_mixture_rows.append(test_sat_state.mixture_with_p(fluid=mixture, sat_line='l', p_units='bar'))
table_mixture_rows.append(test_sat_state.mixture_with_p(fluid=mixture, sat_line='v', p_units='kpa'))
table_mixture_rows.append(test_sat_state.mixture_with_p(fluid=mixture, sat_line='v', p_units='bar'))

# Mixture: output eps_min, eps_max as a pretty table
table_mixture.add_rows(table_mixture_rows)
table_mixture.align["Sat. state at liq/vap curve"] = "l"
table_mixture.align["prop, [units] with : EPSmin, %"] = "c"
table_mixture.align["prop, [units] with : EPSmax, %"] = "c"

print(' ')
print(table_mixture.get_string(title='Saturation state: butane-ethane-methane  '
                                     '/ EPS = [ |prop.model[i] - prop.calc[i]| / prop.model[i]*100 ], % /'))

# My comments on sat_states calculations:
print(' ')
print('MY COMMENTS ON SAT_STATE TESTING: both butane and mixture showed exact fit of model and calculated data. ')
print('for all considered parametersr relative errors made up 0.000000. in the situation when all errors are 0,')
print('the first parameter in the list "t, [k]" was put into the results output table.')
print(' ')
# STOP SAT_STATES TESTING ----------------------------------------------------------------------------------------------

# START SPEC_STATES TESTING --------------------------------------------------------------------------------------------
test_spec_state_tp(mixture)
test_spec_state_td(mixture)
test_spec_state_ph(mixture)
test_spec_state_ps(mixture)
test_spec_state_pd(mixture)
test_spec_state_ds(mixture)
# [HSFLSH error 260] h-s inputs are two-phase or out of bounds, iterative routine is not available to
# find a solution. Test calc. failed in 2ph envelope?!!
test_spec_state_hs(mixture)
test_spec_state_tq(mixture)
# STOP SPEC_STATES TESTING ---------------------------------------------------------------------------------------------

