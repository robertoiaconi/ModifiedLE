# Bisection method for finding m_c such that, at the cut radius, mMLE + m_c = m_h
def bisection(lower, upper, tol, func, *args, **kwargs):
    middle = 0.5 * (lower + upper) # Middle value of mc_on_mh
    while 0.5 * (upper - lower) > tol:

        vals = func(lower, upper, middle, *args, **kwargs)
        lower_criterion = vals[0]
        middle_criterion = vals[1]
        other = vals[2:]

        if middle_criterion == 0: # We're exactly right!
            return (middle, *other)
        elif lower_criterion * middle_criterion < 0: # If the low and middle criteria are of differing signs (i.e. they bracket the zero value)
            upper = middle                                     # then bring the top mc_on_mh down to the middle value
        else:                                                  # else
            lower = middle                                     # bring the bottom value up to the top value
        middle = 0.5 * (lower + upper)
    return (middle, *other)
