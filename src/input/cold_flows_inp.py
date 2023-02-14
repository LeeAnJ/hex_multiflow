"""
этом файле вводим данные горячих потоков: флюид, массовый расход, тем-ра на входе и давление в виде словарей.
словари запаковываем в список flow_in_hot_pack; список flow_in_hot_pack сортируем по убыванию тем-ры потов на входе
"""
from src.rp10.fluid.fluid_class import RP10Fluid
from src.isobar.my_data_classes import Temperature, Pressure, MassRate

# если в потоках один и тот же fluid, то удобнее задавать его отдельно
fluid_0 = RP10Fluid(names=("isobutane", "ethane", "methane"), composition=((0.60, 0.10, 0.30), 'kg/kg'))
# fluid_1 = RP10Fluid(names=("isobutane", "ethane"), composition=((0.60, 0.40), 'kg/kg'))
# fluid_2 = RP10Fluid(names=("isobutane",), composition=((1.00,), 'kg/kg'))

flow_in_cold_pack = list()

# cold flow 1
# ----------------------------------------------------------------------------------------------------------------------
flow_in_cold_pack.append({'fluid': fluid_0,  # or RP10Fluid(names=("isobutane",...), composition=((0.60,...), 'kg/kg'))
                         'mass_rate': MassRate(value=0.009, units='kgs'), # units available: [kgs] or [mols]
                         't_inlet': Temperature(value=143.15, units='k'),    # units available: [k] or [c]
                         'p': Pressure(value=1.8, units='bar')})          # units available: [kpa] or [bar]
# cold flow 2
# ----------------------------------------------------------------------------------------------------------------------
flow_in_cold_pack.append({'fluid': fluid_0,  # or RP10Fluid(names=("isobutane",...), composition=((0.60,...), 'kg/kg'))
                         'mass_rate': MassRate(value=0.009, units='kgs'), # units available: [kgs] or [mols]
                         't_inlet': Temperature(value=153.15, units='k'),    # units available: [k] or [c]
                         'p': Pressure(value=2.0, units='bar')})          # units available: [kpa] or [bar]
# cold flow 3
# ----------------------------------------------------------------------------------------------------------------------
# flow_in_cold_pack.append({'fluid': fluid_0,  # or RP10Fluid(names=("isobutane",...), composition=((0.60,...), 'kg/kg'))
#                          'mass_rate': MassRate(value=0.006, units='kgs'), # units available: [kgs] or [mols]
#                          't_inlet': Temperature(value=163.15, units='k'),    # units available: [k] or [c]
#                          'p': Pressure(value=1.8, units='bar')})         # units available: [kpa] or [bar]

# cold flow 4
# ----------------------------------------------------------------------------------------------------------------------
# flow_in_cold_pack.append({'fluid': fluid_0,  # or RP10Fluid(names=("isobutane",...), composition=((0.60,...), 'kg/kg'))
#                          'mass_rate': MassRate(value=0.005, units='kgs'), # units available: [kgs] or [mols]
#                          't_inlet': Temperature(value=133.15, units='k'),    # units available: [k] or [c]
#                          'p': Pressure(value=1.8, units='bar')})         # units available: [kpa] or [bar]

# print(flow_in_cold_pack[0])
# print(flow_in_cold_pack[1])
# print(flow_in_cold_pack[2])
# print(flow_in_cold_pack[3])
# print('')

# создадим "сортировочную" ф-цию: передаем в нее эл-т списка (а это словарь) и из этого эл-та возвращаем
# параметр, по которому и будем сортировать весь список: list.sort(key=сорт.ф-ция, reverse=True)
def get_t_inlet_as_sorting_parameter(flow: dict) -> float:
    return flow['t_inlet'].value

# отсортируем список flow_in_cold_pack по убыванию тем-ры потока на входе:
flow_in_cold_pack.sort(key=get_t_inlet_as_sorting_parameter, reverse=True)

# print(flow_in_cold_pack[0])
# print(flow_in_cold_pack[1])
# print(flow_in_cold_pack[2])
# print(flow_in_cold_pack[3])
#
# sys.exit('ok')
