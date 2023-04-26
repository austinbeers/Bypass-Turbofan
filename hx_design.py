from scipy.optimize import newton
from atmosphere import isa
from numpy import linspace, zeros, array, meshgrid, shape
from matplotlib.pyplot import plot, show, xlabel, ylabel, title, legend, gca, contour, contourf, colorbar

def hx(x, *args):
    # Extract args
    alt = args[0]  # m
    mach = args[1]
    # Geometry parameters
    A2 = args[2]  # m^2
    ARn = args[3]  # m^2
    # Pressure drop
    deltap = args[4]  # Pa
    # Heat load
    qdot = args[5]  # W

    An = A2 * ARn
    # Inputs
    T1 = x[0]
    p1 = x[1]
    rho1 = x[2]
    u1 = x[3]
    p2 = x[4]
    rho2 = x[5]
    u2 = x[6]
    pn = x[7]
    un = x[8]
    rhon = x[9]

    # Constants
    R = 287
    cp = 1.0  # [J/kg-K]
    gamma = 1.4

    # Ambient conditions
    pa, pta, Ta, Tta, rhoa = isa(alt, mach)
    ua = mach * (gamma * R * Ta)**0.5

    # Diffuser energy conservation
    r1 = cp * Ta + ua**2 / 2 - cp * T1 - u1**2 / 2
    # Isentropic Diffuser
    r2 = (R * Ta)**gamma / pa**(gamma - 1) - p1 / rho1**gamma
    # Continuity across heat exchanger
    r3 = rho1 * u1 - rho2 * u2
    # Momentum across heat exchanger
    r4 = p1 + rho1 * u1**2 - p2 - rho2 * u2**2 - deltap * 0.5 * rho1 * u1**2
    # Energy across heat exchanger
    r5 = (gamma / (gamma - 1)) * p1 / rho1 + u1**2 / 2 - (gamma / (gamma - 1)) * p2 / rho2 - u2**2 / 2 + qdot / rho2 / u2 / A2
    # Heat exchanger to nozzle continuity
    r6 = rho2 * u2 - rhon * un * An / A2
    # Heat exchanger to nozzle energy
    r7 = (gamma / (gamma - 1)) * p2 / rho2 + u2**2 / 2 - (gamma / (gamma - 1))*pn / pn - un**2 / 2
    # Isentropic nozzle
    r8 = p2 / rho2**gamma - pn / rhon**gamma
    r9 = 0
    r10 = 0

    res = [r1, r2, r3, r4, r5, r6, r7, r8, r9, r10]
    return res


if __name__ == "__main__":

    # Geometry parameters
    A2 = 0.8  # m^2
    ARn = 1.0  # m^2
    # Pressure drop
    deltap_var = linspace(5, 12, 5)  # [Pa] Non-dimensionalized pressure drop
    # Sweep heat load
    qdot_var = linspace(50000, 100000, 6)  # W

    qdot_plot, deltap_plot = meshgrid(qdot_var, deltap_var)
    thrust_plot = zeros(shape(qdot_plot))
    mach_plot = zeros(shape(qdot_plot))

    # Generate initial guess
    alt = 10972  # m
    mach = 0.78
    gamma = 1.4
    R = 287
    pa, pta, Ta, Tta, rhoa = isa(alt, mach)
    ua = mach * (gamma * R * Ta) ** 0.5
    x0 = [Ta, pa, rhoa, ua, pa, rhoa, ua, pa, ua, rhoa]
    An = A2 * ARn


    for dp_idx in range(0, len(deltap_var)):
        deltap = deltap_var[dp_idx]
        for q_idx in range(0, len(qdot_var)):
            qdot = qdot_var[q_idx]
            root = newton(hx, x0, args=(alt, mach, A2, ARn, deltap, qdot), tol=1E-8, maxiter=50)
            # Calculate thrust
            rhon = root[9]
            un = root[8]
            pn = root[7]
            Tn = pn / rhon / R
            mdot = rhon * un * An
            thrust = mdot * (un - ua) + A2 * ARn * (pn - pa)
            print('{}, {}'.format(qdot, thrust))
            thrust_plot[dp_idx][q_idx] = thrust / 1000  # [kN]
            mach_plot[dp_idx][q_idx] = un / (gamma * R * Tn)**0.5  # [kN]

    # Plot results contour
    contourf(qdot_var, deltap_var, thrust_plot)
    colorbar()
    # current_values = gca().get_yticks()
    # gca().set_yticklabels(['{:,.5f}'.format(x) for x in current_values])
    title('Thrust')
    xlabel('Heat Load [W]')
    ylabel(r'$\Psi$')
    show()

    # Plot results contour
    contourf(qdot_var, deltap_var, mach_plot)
    colorbar()
    title(r'$M_e$')
    xlabel('Heat Load [W]')
    ylabel(r'$\Psi$')
    show()