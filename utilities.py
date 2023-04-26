from numpy import linspace, interp, array, zeros, shape

def calc_D(mach, gamma):
    # Calculate D(M) for corrected flow per unit area
    D = mach / (1 + 0.5 * (gamma - 1) * mach**2)**(0.5 * ((gamma + 1) / (gamma - 1)))
    return D

def calc_M(D, gamma):
    # Calculate mach number from corrected flow per unit area
    mach_var = linspace(0, 1.0, 200)
    D_m = zeros(shape(mach_var))
    for idx in range(0, len(mach_var)):
        mach = mach_var[idx]
        D_m[idx] = mach / (1 + 0.5 * (gamma - 1) * mach**2)**(0.5 * (gamma + 1) / (gamma - 1))

    # Interpolate mach from list of values
    res = interp(D, D_m, mach_var)
    return res