from compare_calculated_and_sample_data import compare_calculated_and_sample_data


def calc_pd_and_compare_with_sample(fluid=None,
                                    p=None, p_units='kpa',
                                    d=None, d_units='moll',
                                    x=None, x_units='molmol'):

    # f(p, d, x)
    fluid.calc_spec_state(p=(p, p_units), d=(d, d_units), x=(x, x_units))
    return compare_calculated_and_sample_data(fluid=fluid,
                                       x_name='p', x_units=p_units,
                                       y_name='d', y_units=d_units,
                                       z_name='x', z_units=x_units)
