"""
An example of calling the legacy API of REFPROP
By Ian Bell, NIST, 2018, ian.bell@nist.gov
https://github.com/usnistgov/REFPROP-wrappers/blob/master/wrappers/python/ctypes/examples/test_SETREF.py

temperature                     K
pressure, fugacity              kPa
density                         mol/L
composition                     mole fraction
quality                         mole basis (moles vapor/total moles)
enthalpy, internal energy       J/mol
Gibbs, Helmholtz free energy    J/mol
entropy, heat capacity          J/(mol.K)
speed of sound                  m/s
Joule-Thomson coefficient       K/kPa
d(p)/d(rho)                     kPa.L/mol
d2(p)/d(rho)2                   kPa.(L/mol)^2
viscosity                       microPa.s (10^-6 Pa.s)
thermal conductivity            W/(m.K)
dipole moment                   debye
surface tension                 N/m

TRANSPORT PROPERTIES:
   viscosity [microPa.s (10^-6 Pa.s)] and thermal conductivity [W/(m.K)] are calculated for single phase only:
               vis, therm_cond, ierr, herr = r.TRNPRPdll(T,D,z)
   for 2-phase: 0 < q < 1 one should calculate for each phase independently
               vis_liq, therm_cond_liq, ierr, herr = r.TRNPRPdll(T,Dl,x), where x - liq.phase composition
               vis_vap, therm_cond_vap, ierr, herr = r.TRNPRPdll(T,Dv,y), where y - vap.phase composition

   surface tension [N/m] for liquid phase only:
               sigma, ierr, herr = r.SURFTdll(T, Dl, x), where x - liq.phase composition

   Joule-Thomson coefficient [K/kPa]:
               P,e,h,s,Cv,Cp,w,hjt = r.THERMdll(T,D,z)
"""

from __future__ import print_function
from src.rp10.fluid.fluid_class import RP10Fluid

# input fluid's data
fluid_component_name = tuple(["butane"])
# component_composition_mol = tuple([1])
# component_composition_kg = tuple([1])

# mixture = RP10Fluid(("butane",))
# mixture = RP10Fluid(names=("butane",))
# mixture = RP10Fluid(("butane",), composition=((1.0,),))
# mixture = RP10Fluid(names=("butane",), composition=((1.0,), 'mol/mol'))
# mixture = RP10Fluid(names=("butane",), composition=(arr.array('d', [1] + [0]*19), 'mol/mol'))

# mixture = RP10Fluid(("butane",))
mixture = RP10Fluid(names=("butane", "ethane", "methane"), composition=((0.6, 0.3, 0.1), 'mol/mol'))
# mixture = RP10Fluid(names=("butane", "ethane", "methane"), composition=((0.6, 0.3, 0.1), 'mol/mol'))
# mixture = RP10Fluid(names=("butane", "ethane", "methane"), composition=((0.7665, 0.1983, 0.0353), 'kg/kg'))

# mixture.calc_sat_state(sat_curve_flag='v', t=(25.0, 'c'), x=(arr.array('d', [0.4, 0.4, 0.2] + [0]*17), 'kg/kg'))                #, x=((0.4, 0.4, 0.2,), 'molmol'))
# mixture.calc_sat_state(sat_curve_flag='v', t=(25.0, 'c'), x=((0.4, 0.4, 0.2,), 'kg/kg'))                #, x=((0.4, 0.4, 0.2,), 'molmol'))
# mixture.calc_sat_state(sat_curve_flag='v', t=(300.0, 'K'))                #, x=((0.4, 0.4, 0.2,), 'molmol'))
mixture.calc_sat_state(sat_curve_flag='v', p=(100.0, 'kPa'))                #, x=((0.4, 0.4, 0.2,), 'molmol'))

if mixture.error.index > 0:
    mixture.error.print_and_terminate()
else:
    mixture.convert_dataset_to_user_units(flag_data='dew')
    print(mixture.state.data['dew'])
    # mixture.print_spec_state(units_tag='internal')
    # mixture.print_spec_state(units_tag='user')

    # mixture.print_sat_state(sat_curve_symbol='l', units_tag='internal')
    # mixture.print_sat_state(sat_curve_symbol='l', units_tag='user')
