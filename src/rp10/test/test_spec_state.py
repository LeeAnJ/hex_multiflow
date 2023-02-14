from src.rp10.fluid.fluid_class import RP10Fluid
import data.butane_ethane_methane_spec_state as sample
from test_spec_state_hs import test_spec_state_hs
from test_spec_state_pd import test_spec_state_pd
from test_spec_state_ph import test_spec_state_ph
from test_spec_state_ps import test_spec_state_ps
from test_spec_state_td import test_spec_state_td
from test_spec_state_tp import test_spec_state_tp
from test_spec_state_ds import test_spec_state_ds
from test_spec_state_tq import test_spec_state_tq

# initiate RP10Fluid mixture to be tested: ("butane", "ethane", "methane"), (0.6, 0.3, 0.1), [mol/mol].
# get mixture composition from sample dict: sample.rp10_spec_state_liq

x_units = 'molmol'
x = sample.rp10_spec_state_liq['x'][x_units]

mixture = RP10Fluid(names=("butane", "ethane", "methane"), composition=(x, x_units))

# calculate specific state
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
