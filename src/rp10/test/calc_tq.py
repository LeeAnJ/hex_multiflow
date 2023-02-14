from compare_calculated_and_sample_data import compare_calculated_and_sample_data


def calc_tq_and_compare_with_sample(fluid=None,
                                    t=None, t_units='k',
                                    q=None, q_units='molmol',
                                    x=None, x_units='molmol'):

    # f(t, q, x)
    fluid.calc_spec_state(t=(t, t_units), q=(q, q_units), x=(x, x_units))
    return compare_calculated_and_sample_data(fluid=fluid,
                                       x_name='t', x_units=t_units,
                                       y_name='q', y_units=q_units,
                                       z_name='x', z_units=x_units)
