from compare_calculated_and_sample_data import compare_calculated_and_sample_data


def calc_td_and_compare_with_sample(fluid=None,
                                    t=None, t_units='k',
                                    d=None, d_units='moll',
                                    x=None, x_units='molmol'):

    # f(t, d, x)
    fluid.calc_spec_state(t=(t, t_units), d=(d, d_units), x=(x, x_units))
    return compare_calculated_and_sample_data(fluid=fluid,
                                       x_name='t', x_units=t_units,
                                       y_name='d', y_units=d_units,
                                       z_name='x', z_units=x_units)
