from compare_calculated_and_sample_data import compare_calculated_and_sample_data


def calc_hs_and_compare_with_sample(fluid=None,
                                    h=None, h_units='jmol',
                                    s=None, s_units='jmolk',
                                    x=None, x_units='molmol'):
    # f(h, s, x)
    fluid.calc_spec_state(h=(h, h_units), s=(s, s_units), x=(x, x_units))
    return compare_calculated_and_sample_data(fluid=fluid,
                                              x_name='h', x_units=h_units,
                                              y_name='s', y_units=s_units,
                                              z_name='x', z_units=x_units)
