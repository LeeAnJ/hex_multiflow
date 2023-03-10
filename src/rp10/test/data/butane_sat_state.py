import array as arr

# butane = RP10Fluid(names=("butane",))

# satt = f(T=300 K)
satt_liq_model_pattern = {
    't': {'k': 300.0, 'c': 26.850000000000023},
    'p': {'kpa': 257.59575618841797, 'bar': 2.5759575618841795},
    'd': {'moll': 9.818612791950082, 'kgm3': 570.6793764162811},
    'h': {'jmol': 15343.993677193906, 'jkg': 263995.40411742684},
    'e': {'jmol': 15317.758224014653, 'jkg': 263544.0197379771},
    's': {'jmolk': 71.05170291258403, 'jkgk': 1222.4537769145702},
    'cp': {'jmolk': 142.47004096924923, 'jkgk': 2451.2155591021883},
    'cv': {'jmolk': 100.51108393726729, 'jkgk': 1729.3062536735927},
    'q': {'molmol': 0, 'kgkg': 0},
    'w': {'ms': 890.8783373854574},
    'eta': {'upas': 156.97806888522584},
    'tcx': {'wmk': 0.10393141384483205},
    'mm': {'gmol': 58.1222},
    'x': {'molmol': arr.array('d', [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
          'kgkg': arr.array('d', [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])}
}
satt_vap_model_pattern = {
    't': {'k': 300.0, 'c': 26.850000000000023},
    'p': {'kpa': 257.59575618841797, 'bar': 2.5759575618841795},
    'd': {'moll': 0.1121152347831493, 'kgm3': 6.516384099113161},
    'h': {'jmol': 36243.57446804926, 'jkg': 623575.406093528},
    'e': {'jmol': 33945.97631658663, 'jkg': 584044.9314820607},
    's': {'jmolk': 140.71697221543522, 'jkgk': 2421.053783501575},
    'cp': {'jmolk': 105.26664124939859, 'jkgk': 1811.12623488785},
    'cv': {'jmolk': 93.09826727907578, 'jkgk': 1601.7677802814724},
    'q': {'molmol': 1, 'kgkg': 1},
    'w': {'ms': 202.15069202023003},
    'eta': {'upas': 7.377097828275173},
    'tcx': {'wmk': 0.016784317058438067},
    'mm': {'gmol': 58.1222},
    'x': {'molmol': arr.array('d', [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
          'kgkg': arr.array('d', [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])}}

# satp = f(p=100 kPa)
satp_liq_model_pattern = {
    't': {'k': 272.31397194261064, 'c': -0.8360280573893419},
    'p': {'kpa': 100.0, 'bar': 1.0},
    'd': {'moll': 10.351114599463703, 'kgm3': 601.6295529729492},
    'h': {'jmol': 11512.037501328734, 'jkg': 198066.10041135288},
    'e': {'jmol': 11502.376705958393, 'jkg': 197899.88517224733},
    's': {'jmolk': 57.71120793018112, 'jkgk': 992.9288280584891},
    'cp': {'jmolk': 134.14501609685033, 'jkgk': 2307.9824249056355},
    'cv': {'jmolk': 95.15741876507946, 'jkgk': 1637.1957490439015},
    'q': {'molmol': 0, 'kgkg': 0},
    'w': {'ms': 1039.7008327160645},
    'eta': {'upas': 204.10289982962763},
    'tcx': {'wmk': 0.11567876589562054},
    'mm': {'gmol': 58.1222},
    'x': {'molmol': arr.array('d', [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
          'kgkg': arr.array('d', [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])}
}

satp_vap_model_pattern = {
    't': {'k': 272.31397194261064, 'c': -0.8360280573893419},
    'p': {'kpa': 100.0, 'bar': 1.0},
    'd': {'moll': 0.0460447386201796, 'kgm3': 2.6762215070298025},
    'h': {'jmol': 33948.15487865798, 'jkg': 584082.4139254532},
    'e': {'jmol': 31776.354082373866, 'jkg': 546716.2991485846},
    's': {'jmolk': 140.10182937819698, 'jkgk': 2410.470171091201},
    'cp': {'jmolk': 95.26691718837289, 'jkgk': 1639.0796836384873},
    'cv': {'jmolk': 85.13136392718064, 'jkgk': 1464.6961733585558},
    'q': {'molmol': 1, 'kgkg': 1},
    'w': {'ms': 200.0743056246413},
    'eta': {'upas': 6.718666571853523},
    'tcx': {'wmk': 0.014113489994994463},
    'mm': {'gmol': 58.1222},
    'x': {'molmol': arr.array('d', [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
          'kgkg': arr.array('d', [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])}
}

