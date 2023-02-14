import random
import sys
import timeit

# import cProfile, pstats, io
# from pstats import SortKey

from src.rp10.fluid.fluid_class import RP10Fluid
from src.isobar.isobar import Isobar
from src.isobar.my_data_classes import Pressure, Temperature

# input fluid's data:
mixture = RP10Fluid(names=("isobutane", "ethane", "methane"),
                    composition=((0.60, 0.10, 0.30), 'kg/kg'))

# input isobar parameters:
t_min_k = 153.15
t_max_k = 333.15
p_high = Isobar(fluid=mixture,
                p=Pressure(value=20.0, units='bar'),
                t_min=Temperature(value=t_min_k, units='k'),
                t_max=Temperature(value=t_max_k, units='k'))

# for i in range(len(p_high.t_k)):
#     print(p_high.h_jmol[i], p_high.t_k[i])
# print(len(p_high.t_k))
# sys.exit()

# создадим списки "базовых" t[i] и h[i] точек (они же будут аргументами расчет. функция) на изобаре
# случайным образом в заданных пределах и общим числом n_rand:
n_rand = 1000
t_random = [random.uniform(t_min_k, t_max_k) for i in range(n_rand)]
h_random = [random.uniform(p_high.h_min_jmol, p_high.h_max_jmol) for j in range(n_rand)]
# for i in range(n_rand):
#     print(i, t_random[i], h_random[i])
# sys.exit()

n_tau_iterations = 100     # repeat t_iterations or h_iterations n_tau_iterations times
n_t_iterations = 600       # number of temperatures for testing within n_rand: n_t_iterations <= n_rand
n_h_iterations = 600       # number of enthalpies for testing within n_rand: n_h_iterations <= n_rand

# init "errors" variables
err_t_min = float('+inf')
err_t_max = float('-inf')
err_t_avr = 0
err_t_sum = 0

err_s_min = float('+inf')
err_s_max = float('-inf')
err_s_avr = 0
err_s_sum = 0

# pr = cProfile.Profile()
# pr.enable()

# for j in range(n_tau_iterations):     # это верх. цикл, если я хочу провести усреднение еще и по времени
start = timeit.default_timer()

# for i in range(n_t_iterations):
for i in range(n_h_iterations):

    t_linear, s_linear = p_high.get_t_k_s_jmolk_with_h__linear_interpolation(h_jmol=h_random[i])
    t_rp_10, s_rp_10 = p_high.get_t_k_s_jmolk_with_h__rp10(h_jmol=h_random[i])

    err_t = abs(t_linear - t_rp_10)/t_rp_10 * 100.0
    err_s = abs(s_linear - s_rp_10) /s_rp_10 * 100.0

    if err_t <= err_t_min:
        err_t_min = err_t
    if err_t >= err_t_max:
        err_t_max = err_t
    err_t_sum += err_t

    if err_s <= err_s_min:
        err_s_min = err_s
    if err_s >= err_s_max:
        err_s_max = err_s
    err_s_sum += err_s

    # print(f'{i}\t\t{h_random[i]:10.4f}\t{t_linear:7.2f}\t{s_linear:10.4f}\t{err_t:8.5f}\t{err_s:8.5f}')

# err_avr = err_sum/n_t_iterations
err_t_avr = err_t_sum/n_h_iterations
err_s_avr = err_s_sum/n_h_iterations

print('---------------------------------------------------------')
print('Relative error: |t_linear_interp -t_rp10|/t_rp10 * 100')
print('case: t, s = f(h, p=const)')
print('---------------------------------------------------------')
print(f'err_t_min, % = {err_t_min:9.6f}\t\terr_s_min, % = {err_s_min:9.6f}')
print(f'err_t_max, % = {err_t_max:9.6f}\t\terr_s_max, % = {err_s_max:9.6f}')
print(f'err_t_avr, % = {err_t_avr:9.6f}\t\terr_s_avr, % = {err_s_avr:9.6f}')
print('---------------------------------------------------------')

stop = timeit.default_timer()
dtau_interpol = stop-start
print(f'Время расчета {n_h_iterations} циклов сравнения t_linear_interp с t_rp10, сек: {dtau_interpol:9.6f}') #dtau_avr)  # время расчета см. выше
print(f'Cреднее время расчета на один цикл, сек: {dtau_interpol/n_h_iterations:9.6f}\n\n') #dtau_avr)  # время расчета см. выше

# init "errors" variables
err_h_min = float('+inf')
err_h_max = float('-inf')
err_h_avr = 0
err_h_sum = 0

err_s_min = float('+inf')
err_s_max = float('-inf')
err_s_avr = 0
err_s_sum = 0

# pr = cProfile.Profile()
# pr.enable()

# for j in range(n_tau_iterations):     # это верх. цикл, если я хочу провести усреднение еще и по времени
start = timeit.default_timer()

# for i in range(n_t_iterations):
for i in range(n_t_iterations):

    h_linear, s_linear = p_high.get_h_jmol_s_jmolk_with_t__linear_interpolation(t_k=t_random[i])
    h_rp_10, s_rp_10 = p_high.get_h_jmol_s_jmolk_with_t__rp10(t_k=t_random[i])

    err_h = abs(h_linear - h_rp_10)/h_rp_10 * 100.0
    err_s = abs(s_linear - s_rp_10) /s_rp_10 * 100.0

    if err_h <= err_h_min:
        err_h_min = err_h
    if err_h >= err_h_max:
        err_h_max = err_h
    err_h_sum += err_h

    if err_s <= err_s_min:
        err_s_min = err_s
    if err_s >= err_s_max:
        err_s_max = err_s
    err_s_sum += err_s

    # print(f'{i}\t\t{h_random[i]:10.4f}\t{t_linear:7.2f}\t{s_linear:10.4f}\t{err_t:8.5f}\t{err_s:8.5f}')

# err_avr = err_sum/n_t_iterations
err_h_avr = err_h_sum/n_t_iterations
err_s_avr = err_s_sum/n_t_iterations

print('---------------------------------------------------------')
print('Relative error: |h_linear_interp -h_rp10|/h_rp10 * 100')
print('case: h, s = f(t, p=const)')
print('---------------------------------------------------------')
print(f'err_h_min, % = {err_h_min:10.6f}\t\terr_s_min, % = {err_s_min:9.6f}')
print(f'err_h_max, % = {err_h_max:10.6f}\t\terr_s_max, % = {err_s_max:9.6f}')
print(f'err_h_avr, % = {err_h_avr:10.6f}\t\terr_s_avr, % = {err_s_avr:9.6f}')
print('---------------------------------------------------------')

stop = timeit.default_timer()
dtau_interpol = stop-start
print(f'Время расчета {n_t_iterations} циклов сравнения h_linear_interp с h_rp10, сек: {dtau_interpol:9.6f}') #dtau_avr)  # время расчета см. выше
print(f'Cреднее время расчета на один цикл, сек: {dtau_interpol/n_h_iterations:9.6f}') #dtau_avr)  # время расчета см. выше

#         dtau_sum += dtau_interpol
# dtau_avr = dtau_sum/tau_range/t_range

# pr.disable()
# s=io.StringIO()
# sortby=SortKey.CUMULATIVE
# ps=pstats.Stats(pr,stream=s).sort_stats(sortby)
# ps.print_stats()
# print(s.getvalue())

# start = timeit.default_timer()
# for i in range(10):
#     p_high.get_h_jmol_with_refprop10(t_k=t_random[i])
# stop = timeit.default_timer()
# dtau_rp10 = stop - start
# print('\nВремя расчета, сек: ', dtau_rp10)  # время расчета см. выше
# print('dtau_rp10/dtau_interpol = ', dtau_rp10/dtau_interpol)
#
# print('ok')
# print(t_random[0])
# print(t_random[0], p_high.get_h_jmol(t_k=t_random[0]))
