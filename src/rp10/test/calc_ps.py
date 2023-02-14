from compare_calculated_and_sample_data import compare_calculated_and_sample_data


def calc_ps_and_compare_with_sample(fluid=None,
                                    p=None, p_units='kpa',
                                    s=None, s_units='jmolk',
                                    x=None, x_units='molmol'):

    # f(p, s, x)
    fluid.calc_spec_state(p=(p, p_units), s=(s, s_units), x=(x, x_units))
    return compare_calculated_and_sample_data(fluid=fluid,
                                       x_name='p', x_units=p_units,
                                       y_name='s', y_units=s_units,
                                       z_name='x', z_units=x_units)
