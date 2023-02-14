from compare_calculated_and_sample_data import compare_calculated_and_sample_data


def calc_ph_and_compare_with_sample(fluid=None,
                                    p=None, p_units='kpa',
                                    h=None, h_units='jmol',
                                    x=None, x_units='molmol'):

    # f(p, h, x)
    fluid.calc_spec_state(p=(p, p_units), h=(h, h_units), x=(x, x_units))
    return compare_calculated_and_sample_data(fluid=fluid,
                                       x_name='p', x_units=p_units,
                                       y_name='h', y_units=h_units,
                                       z_name='x', z_units=x_units)
