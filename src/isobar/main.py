import random

from src.rp10.fluid.fluid_class import RP10Fluid
from src.isobar.isobar import Isobar
from src.isobar.my_data_classes import Pressure, Temperature


if __name__ == "__main__":

    # input fluid's data:
    mixture = RP10Fluid(names=("isobutane", "ethane", "methane"), composition=((0.60, 0.10, 0.30), 'kg/kg'))

    # для построения изобары требуется:
    #       собственно раб.флюид - экз. класса RP10Fluid
    #       темп.диапазон: ]t_min_k .. t_max_k[ ВНИМАНИЕ: на границах диапазона лин. интерп. не работает
    #       давление изобары: экз. класса Pressure
    t_min_k = 153.15
    t_max_k = 333.15
    # cоздадим экз. класса Isobar
    p_high = Isobar(fluid=mixture,
                    p=Pressure(value=20.0, units='bar'),
                    t_min=Temperature(value=t_min_k, units='k'),
                    t_max=Temperature(value=t_max_k, units='k'))

    # выберем пару "базовых" точек (они же будут аргументами расчет. функция) на изобаре
    t_random = random.uniform(t_min_k, t_max_k)
    h_random = random.uniform(p_high.h_min_jmol, p_high.h_max_jmol)

    print(f'\nrandom t,[K] = {t_random} within t_min_k = {t_min_k} .. t_max_k = {t_max_k}')
    print(f'random h,[J/mol] = {h_random} within h_min_jmol = {p_high.h_min_jmol} .. h_max_jmol = {p_high.h_max_jmol}\n')

    # расчет т/д св-в fluid для экз. класса Isobar проводится путем лин. интерполяции между базовыми точками изобары;
    # независимый агрумент может быть: либо t,[K], либо h, [J/mol]

    # h,s = f(t) along isobar.p_kpa
    h_1, s_1 = p_high.get_h_jmol_s_jmolk_with_t__linear_interpolation(t_k=t_random)
    # s = f(t) along isobar.p_kpa
    s_2 = p_high.get_s_jmolk_with_t__linear_interpolation(t_k=t_random)
    # h = f(t) along isobar.p_kpa
    h_2 = p_high.get_h_jmol_with_t__linear_interpolation(t_k=t_random)
    print(t_random, h_1, s_1)
    print(t_random, h_2, s_2)

    # t,s = f(h) along isobar.p_kpa
    t_3, s_3 = p_high.get_t_k_s_jmolk_with_h__linear_interpolation(h_jmol=h_random)
    # t = f(h) along isobar.p_kpa
    t_4 = p_high.get_t_k_with_h_jmol__linear_interpolation(h_jmol=h_random)
    print(h_random, t_3, s_3)
    print(h_random, t_4)

    # также в классе Isobar присутствуют 2 метода для расчета t_k, h_jmol и s_jmolk непосредственно
    # с помощью ф-ций RefProp 10.0 (вероятно, я их использовал для поверки в процессе тестирования изобары)
    h_5, s_5 = p_high.get_h_jmol_s_jmolk_with_t__rp10(t_k=t_random)
    t_5 = p_high.get_t_k_with_refprop10(h_jmol=h_random)
    print(t_random, h_5, s_5)
    print(t_random, h_5)
