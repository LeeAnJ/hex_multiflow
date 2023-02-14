Расчет ОБРАЗЦОВЫХ (MODEL) свойств чистого вещества и зеотропной смеси проводился в программе RefProp 10.0
для:
		1.	чистого butane
		2.	зеотропной смеси: names = ("butane", "ethane", "methane"), composition = ( (0.6, 0.3, 0.1), 'mol/mol' )

Результаты расчетов копи-пастились непосредственно из табл. граф. интерфейса RefProp 10.0 в  виде словарей в 
py-модули, размещенные в папке src\rp10\test\data:
		sat curve:
			butane_sat_state.py:		satt:  sat_liq = f(T=300 K); sat_vap = f(T=300 K); 
													satp: sat_liq = f(p=100 kPa); sat_vap = f(p=100 kPa) 
			butane_ethane_methane_sat_state.py: 	satt:  sat_liq = f(T=300 K); sat_vap = f(T=300 K); 
																				satp: sat_liq = f(p=100 kPa); sat_vap = f(p=100 kPa) 		
																				
		single phase: liq, vap;	2-phase: 2ph_liq, 2ph_blk, 2ph_vap	{blk stays for "bulk"}
			только для смеси (уже и не помню почему для бутана не сделал)
			butane_ethane_methane_spec_state.py: 
																					1ph_liq = f(t=120, K; p=100, kPa)
																						2ph_blk = f(t=200, K; p=100, kPa)
																						2ph_liq = f(t=200, K; p=100, kPa)
																						2ph_vap = f(t=200, K; p=100, kPa)
																							1ph_vap = f(t=270, K; p=100, kPa)
