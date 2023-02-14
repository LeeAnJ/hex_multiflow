from compare_calculated_and_sample_data import compare_calculated_and_sample_data


def calc_tp_and_compare_with_sample(fluid=None,
                                    t=None, t_units='k',
                                    p=None, p_units='kpa',
                                    x=None, x_units='molmol'):

    # f(t, p, x)
    fluid.calc_spec_state(t=(t, t_units), p=(p, p_units), x=(x, x_units))
    return compare_calculated_and_sample_data(fluid=fluid,
                                       x_name='t', x_units=t_units,
                                       y_name='p', y_units=p_units,
                                       z_name='x', z_units=x_units)
