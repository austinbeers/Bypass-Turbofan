from scipy.optimize import newton, root, root_scalar
from atmosphere import isa
from math import e
from utilities import calc_M, calc_D
from numpy import linspace, zeros, array, meshgrid, shape
from matplotlib.pyplot import plot, show, xlabel, ylabel, title, legend, gca, contour, contourf, colorbar

def hx(x, *args):
    # Extract args
    alt = args[0]  # m
    mach = args[1]
    # Geometry parameters
    A = args[2]  # m^2
    beta = args[3]  # [1/m]
    m_dot_h = args[4]

    # Inputs
    mdot = x[0]
    T1 = x[1]
    p1 = x[2]
    u1 = x[3]
    rho1 = x[4]
    T2 = x[5]
    p2 = x[6]
    u2 = x[7]
    rho2 = x[8]
    T02 = x[9]
    p02 = x[10]
    Te = x[11]
    pe = x[12]
    ue = x[13]
    rhoe = x[14]
    NTU = x[15]
    epsilon = x[16]
    Q = x[17]


    # Constants
    R = 287  # [J/kg-K]
    cp = 1004  # [J/kg-K]
    gamma = 1.4

    # Ambient conditions
    pa, p0a, Ta, T0a, rho = isa(alt, mach)

    T01 = T0a
    p01 = p0a
    # Unknowns
    # General - mdot
    # Diffuser - T1, P1, u1, rho1
    # HX -  T2, P2, u2, rho2, T02, P02,
    # Nozzle - Te, Pe, ue, rhoe

    # Isentropic diffuser
    r0 = T01 - T1 - u1**2 / 2 / cp
    r1 = p01 - p1 - 0.5 * rho1 * u1**2
    r2 = p1 - rho1 * R * T1
    r3 = mdot - rho1 * u1 * A

    # Heat exchanger design
    # Heat Addition
    Uh = 385  # [] Conductivity of copper
    Uc = 0.024  # [] Conductivity of air - should this be copper as well?
    Uavg = (Uh + Uc) / 2  # Average Conductivity
    # beta = 500  # [1/m]

    Cph = 2200  # [J/kg-K] Specific heat of glycol
    Cpc = cp  # [J/kg-K] Specific heat of air

    Tc_in = T1  # [K] Inlet temperature
    Th_in = 394  # [K] Temperature of working fluid

    Ch = m_dot_h * Cph

    Cc = mdot * Cpc

    if Cc > Ch:
        Cmax = Cc
        Cmin = Ch
    else:
        Cmax = Ch
        Cmin = Cc

    # Cmin_Cmax = Cmin / Cmax
    r15 = NTU - A * Uavg * beta / Cmin
    r16 = epsilon - (1 - e**(-1*NTU * (1 - Cmin / Cmax)))/(1 - (Cmin / Cmax) * e**(-1 * NTU * (1 - Cmin / Cmax)))
    Qmax = Cmin * (Th_in - Tc_in)
    r17 = Q - epsilon * Qmax

    # Pressure Drop
    sigma = 0.9  # notional guess
    mu1 = 1.458E-6 * T1 ** 1.5 / (T1 + 110.4)  # [Pa*s] Dynamic viscosity at station 2
    mu2 = 1.458E-6 * T2 ** 1.5 / (T2 + 110.4)  # [Pa*s] Dynamic viscosity at station 1
    # Calculate rh from bete
    rh = 0.0254 * 33.7 * beta ** (-0.993)
    Pr = 0.72
    # Calculate l_rh
    l_rh = 0.8 ** 2 * (Pr * rho1 * u1 * rh) / (sigma * mu1)
    l = l_rh * rh  # [m]
    Re = rho1 * u1 * l / mu1
    Cf = 1.328 / Re ** 0.5
    Pf = Cf * l_rh / sigma ** 2
    phi = (1 + mu2 * u2 / mu1 / u1) * Pf / 2 + 2 * (u2 / u1 - 1)

    r4 = T2 - T1 - Q / mdot / Cpc
    r5 = T02 - T2 - u2**2 / cp
    r6 = mdot - rho2 * u2 * A # Continuity through screen
    r7 = p2 - p1 + 0.5 * rho1 * u1 ** 2 * phi
    r8 = p02 - p2 - 0.5 * rho2 * u2**2
    r9 = p2 - rho2 * R * T2

    # Isentropic nozzle
    T0e = T02
    p0e = p02
    r10 = pe - pa
    Me = ((2 / (gamma - 1)) * ((p0e / pe)**((gamma - 1) / gamma) - 1))**0.5
    if Me > 1:  # check if nozzle is choked
        print('Fan nozzle is choked.')
        Me = 1
        r10 = pe - p0e * (1 + ((gamma - 1) / 2) * Me**2)**(-1 * gamma / (gamma - 1))
    else:
        print('Fan nozzle is not choked. M_e = {}'.format(Me))

    r11 = Te - T0e / (1 + ((gamma - 1) / 2) * Me**2)
    r12 = ue - Me * (gamma * R * Te)**0.5
    r13 = mdot - rhoe * ue * A
    r14 = pe - rhoe * R * Te

    # HX Parameters - r15, r16, r17

    res = [r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17]
    return res


def hx_sys(x, *args):
    # Extract args
    alt = args[0]  # m
    mach = args[1]
    # Geometry parameters
    A = args[2]  # m^2
    beta = args[3]  # [1/m]
    Q_sys = args[4]

    # Inputs
    mdot = x[0]
    T1 = x[1]
    p1 = x[2]
    u1 = x[3]
    rho1 = x[4]
    T2 = x[5]
    p2 = x[6]
    u2 = x[7]
    rho2 = x[8]
    T02 = x[9]
    p02 = x[10]
    Te = x[11]
    pe = x[12]
    ue = x[13]
    rhoe = x[14]
    NTU = x[15]
    epsilon = x[16]
    Th_in = x[17]


    # Constants
    R = 287  # [J/kg-K]
    cp = 1004  # [J/kg-K]
    gamma = 1.4

    # Ambient conditions
    pa, p0a, Ta, T0a, rho = isa(alt, mach)

    T01 = T0a
    p01 = p0a
    # Unknowns
    # General - mdot
    # Diffuser - T1, P1, u1, rho1
    # HX -  T2, P2, u2, rho2, T02, P02,
    # Nozzle - Te, Pe, ue, rhoe

    # Isentropic diffuser
    r0 = T01 - T1 - u1**2 / 2 / cp
    r1 = p01 - p1 - 0.5 * rho1 * u1**2
    r2 = p1 - rho1 * R * T1
    r3 = mdot - rho1 * u1 * A

    # Heat exchanger design
    # Heat Addition
    Uh = 385  # [] Conductivity of copper
    Uc = 0.024  # [] Conductivity of air - should this be copper as well?
    Uavg = (Uh + Uc) / 2  # Average Conductivity
    # beta = 500  # [1/m]

    Cph = 2200  # [J/kg-K] Specific heat of glycol
    Cpc = cp  # [J/kg-K] Specific heat of air

    Tc_in = T1  # [K] Inlet temperature
    # Th_in = 394  # [K] Temperature of working fluid

    m_dot_h = 1  # [kg/s] Mass flow rate of working fluid
    Ch = m_dot_h * Cph

    Cc = mdot * Cpc

    if Cc > Ch:
        Cmax = Cc
        Cmin = Ch
    else:
        Cmax = Ch
        Cmin = Cc


    # Cmin_Cmax = Cmin / Cmax
    r15 = NTU - A * Uavg * beta / Cmin
    r16 = epsilon - (1 - e**(-1*NTU * (1 - Cmin / Cmax)))/(1 - (Cmin / Cmax) * e**(-1 * NTU * (1 - Cmin / Cmax)))
    r17 = Q_sys - epsilon * Cmin * (Th_in - Tc_in)


    # Pressure Drop
    sigma = 0.9  # notional guess
    mu1 = 1.458E-6 * T1 ** 1.5 / (T1 + 110.4)  # [Pa*s] Dynamic viscosity at station 2
    mu2 = 1.458E-6 * T2 ** 1.5 / (T2 + 110.4)  # [Pa*s] Dynamic viscosity at station 1
    # Calculate rh from bete
    rh = 0.0254 * 33.7 * beta**(-0.993)
    Pr = 0.72
    # Calculate l_rh
    l_rh = 0.8**2 * (Pr * rho1 * u1 * rh) / (sigma * mu1)
    l = l_rh * rh  # [m]

    Re = rho1 * u1 * l / mu1
    Cf = 1.328 / Re ** 0.5
    Pf = Cf * l_rh / sigma ** 2
    phi = (1 + mu2 * u2 / mu1 / u1) * Pf / 2 + 2 * (u2 / u1 - 1)

    r4 = T2 - T1 - Q_sys / mdot / Cpc
    r5 = T02 - T2 - u2**2 / cp
    r6 = mdot - rho2 * u2 * A # Continuity through screen
    r7 = p2 - p1 + 0.5 * rho1 * u1 ** 2 * phi
    r8 = p02 - p2 - 0.5 * rho2 * u2**2
    r9 = p2 - rho2 * R * T2

    # Isentropic nozzle
    T0e = T02
    p0e = p02
    r10 = pe - pa
    Me = ((2 / (gamma - 1)) * ((p0e / pe)**((gamma - 1) / gamma) - 1))**0.5
    if Me > 1:  # check if nozzle is choked
        print('Fan nozzle is choked.')
        Me = 1
        r10 = pe - p0e * (1 + ((gamma - 1) / 2) * Me**2)**(-1 * gamma / (gamma - 1))
    else:
        print('Fan nozzle is not choked. M_e = {}'.format(Me))

    r11 = Te - T0e / (1 + ((gamma - 1) / 2) * Me**2)
    r12 = ue - Me * (gamma * R * Te)**0.5
    r13 = mdot - rhoe * ue * A
    r14 = pe - rhoe * R * Te

    # HX Parameters - r15, r16, r17

    res = [r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17]
    return res


if __name__ == "__main__":

    # Generate initial guess
    alt = 10972  # m
    mach = 0.78
    gamma = 1.4
    cp = 1004  # [J/kg-K]
    R = 287
    pa, p0a, Ta, T0a, rhoa = isa(alt, mach)
    ua = mach * (gamma * R * Ta) ** 0.5

    beta_var = linspace(100, 800, 40)  # [1/m^2]
    area_var = linspace(0.05, 0.3, 50)  # [m^2]

    area_plot, beta_plot = meshgrid(area_var, beta_var)
    epsilon_plot = zeros(shape(area_plot))
    q_plot = zeros(shape(area_plot))
    mache_plot = zeros(shape(area_plot))
    phi_plot = zeros(shape(area_plot))
    m_dot_h = 2  # [kg/s]
    for A_idx in range(0, len(area_var)):
        A_hx = area_var[A_idx]  # [m^2]
        for beta_idx in range(0, len(beta_var)):
            beta = beta_var[beta_idx]
            # Inputs -  [mdot, T1, p1, u1, rho1, T2, p2, u2, rho2, T02, p02, Te, pe, ue, rhoe, NTU, epsilon, Q]
            x0 = [rhoa * ua * A_hx, Ta, pa, ua, rhoa, 1.2*Ta, 0.8*pa, 1.2*ua, 0.8*rhoa, 1.2*T0a, 0.8*p0a, Ta, pa, ua, rhoa, 4.0, 0.9, 100000]
            sol = root(hx, x0, args=(alt, mach, A_hx, beta, m_dot_h), tol=1E-8)

            # Extract Outputs
            mdot = sol.x[0]
            T1 = sol.x[1]
            p1 = sol.x[2]
            u1 = sol.x[3]
            rho1 = sol.x[4]
            T2 = sol.x[5]
            p2 = sol.x[6]
            u2 = sol.x[7]
            rho2 = sol.x[8]
            T02 = sol.x[9]
            p02 = sol.x[10]
            Te = sol.x[11]
            pe = sol.x[12]
            ue = sol.x[13]
            rhoe = sol.x[14]
            NTU = sol.x[15]
            epsilon = sol.x[16]
            q = sol.x[17]

            # Net Force
            f_net = mdot * (ue - ua) + A_hx * (pe - pa)
            print('Net Force: {} N'.format(f_net))

            # Mach at Exit
            D = mdot * (R * T02)**0.5 / A_hx / p02 / gamma**0.5
            M_e = calc_M(D, gamma)
            print('Exit Mach: {}'.format(M_e))

            q_plot[beta_idx][A_idx] = q / 1000  # [kW]
            mache_plot[beta_idx][A_idx] = M_e  # [-]
            phi_plot[beta_idx][A_idx] = (p1 - p2) / (0.5 * rho1 * u1**2)

    # Plot Contours
    contourf(area_plot, beta_plot, mache_plot)
    colorbar()
    xlabel(r'Area [$m^2$]')
    ylabel(r'$\beta$ [1/$m^2$]')
    title(r'Exit Mach, $M_e$, for $\dot{m}_h = 2 kg/s$')
    show()

    contourf(area_plot, beta_plot, phi_plot)
    colorbar()
    xlabel(r'Area [$m^2$]')
    ylabel(r'$\beta$ [1/$m$]')
    title(r'Pressure Drop, $\Phi$, for $\dot{m}_h = 2 kg/s$')
    show()

