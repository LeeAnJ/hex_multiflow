from compare_calculated_and_sample_data import compare_calculated_and_sample_data


def calc_ds_and_compare_with_sample(fluid=None,
                                    d=None, d_units='moll',
                                    s=None, s_units='jmolk',
                                    x=None, x_units='molmol'):

    # f(d, s, x)
    fluid.calc_spec_state(d=(d, d_units), s=(s, s_units), x=(x, x_units))
    return compare_calculated_and_sample_data(fluid=fluid,
                                       x_name='d', x_units=d_units,
                                       y_name='s', y_units=s_units,
                                       z_name='x', z_units=x_units)
