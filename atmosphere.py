from math import exp
from numpy import linspace, array, concatenate, tan, append, flip, sum, multiply, mean
from matplotlib.pyplot import plot, show, axis, legend, figure, savefig, subplots, xlabel, ylabel

# inputs:
# h: Altitude in meters
# M: Mach number
def isa(h, M):

    # International Standard Atmosphere
    # (dp/dh) = -(rho)*g
    # p = (rho)*R*T
    #
    # dp/p = -g*dh/(R*T)
    # T = a+b*h
    #
    # b != 0:
    # ln(p2/p1) = -g/(b*R)*ln((a+b*h2)/(a+b*h1))
    # p2 = p1*((a+b*h2)/(a+b*h1))^(-g/(b*R))
    #
    # b == 0:
    # ln(p2/p1) = -g*(h2-h1)/R/T
    # p2 = p1*exp(-g*(h2-h1)/R/T)
    
    g = 9.80665
    R = 287.058
    gamma = 1.4
    T_sls = 288.15  # K

    if h <= 11000:
        # Troposphere
        # a = 19 + 273  # Base temperature (K)
        # b = -0.0065  # Temperature lapse rate (K/m)
        # p0 = 108900  # Base pressure (Pa)
        # h0 = -610  # Base altitude (m)
        a = 288.16  # Base temperature (K)
        b = -0.0065  # Temperature lapse rate (K/m)
        p0 = 101325  # Base pressure (Pa)
        h0 = 0  # Base altitude (m)
        p = p0 * ((a + b * h) / (a + b * h0)) ** (-g / b / R)
        # Temperature
        h_ref = 0
        T_ref = T_sls
        dTdh = -6.5 / 1000  # K/m
        T = T_ref + dTdh * (h - h_ref)

    elif h <= 20000:
        # Tropopause
        a = -56.5 + 273  # Base temperature
        b = 0  # Temperature lapse rate
        p0 = 22632  # Base pressure
        h0 = 11000  # Base altitude
        p = p0 * exp(-g*(h-h0)/R/a)
        p_ref, p_t_ref, T_ref, T_t_ref, rho = isa(h0, 0)
        T = T_ref
    else:
        # Stratosphere (up to 32000 m)
        a = -56.5 + 273 # Base temperature
        b = 0.001 # Temperature lapse rate
        p0 = 5474.9 # Base pressure
        h0 = 20000 # Base altitude
        p = p0 * ((a+b*h)/(a+b*h0))**(-g/b/R)
        p_ref, p_t_ref, T_ref, T_t_ref, rho = isa(h0, 0)
        dTdh = 1 / 1000  # K/m
        T = T_ref + dTdh * (h - h0)
        p = p_ref * (T / T_ref) ** (-g / R / dTdh)

    rho = p / R / T
    # T = a + b * (h-h0)
    p_t = p*(1+(gamma-1)/2*M**2)**(gamma/(gamma-1))
    T_t = T*(1+(gamma-1)/2*M**2)
    
    return p, p_t, T, T_t, rho


def tropopause(h, M):
    # International Standard Atmosphere
    # (dp/dh) = -(rho)*g
    # p = (rho)*R*T
    #
    # dp/p = -g*dh/(R*T)
    # T = a+b*h
    #
    # b != 0:
    # ln(p2/p1) = -g/(b*R)*ln((a+b*h2)/(a+b*h1))
    # p2 = p1*((a+b*h2)/(a+b*h1))^(-g/(b*R))
    #
    # b == 0:
    # ln(p2/p1) = -g*(h2-h1)/R/T
    # p2 = p1*exp(-g*(h2-h1)/R/T)

    g = 9.80665
    R = 287.058
    gamma = 1.4
    T_sls = 288.15  # K

    # Troposphere
    a = -56.5 + 273  # Base temperature
    b = 0  # Temperature lapse rate
    p0 = 22632  # Base pressure
    h0 = 11000  # Base altitude
    p = p0 * exp(-g * (h - h0) / R / a)
    p_ref, p_t_ref, T_ref, T_t_ref, rho = isa(h0, 0)
    T = T_ref

    rho = p / R / T
    T = a + b * (h - h0)
    p_t = p * (1 + (gamma - 1) / 2 * M ** 2) ** (gamma / (gamma - 1))
    T_t = T * (1 + (gamma - 1) / 2 * M ** 2)

    return p, p_t, T, T_t, rho

def troposphere(h, M):
    # International Standard Atmosphere
    # (dp/dh) = -(rho)*g
    # p = (rho)*R*T
    #
    # dp/p = -g*dh/(R*T)
    # T = a+b*h
    #
    # b != 0:
    # ln(p2/p1) = -g/(b*R)*ln((a+b*h2)/(a+b*h1))
    # p2 = p1*((a+b*h2)/(a+b*h1))^(-g/(b*R))
    #
    # b == 0:
    # ln(p2/p1) = -g*(h2-h1)/R/T
    # p2 = p1*exp(-g*(h2-h1)/R/T)

    g = 9.80665
    R = 287.058
    gamma = 1.4
    T_sls = 288.15  # K

    # Troposphere
    a = 288.16  # Base temperature (K)
    b = -0.0065  # Temperature lapse rate (K/m)
    p0 = 101325  # Base pressure (Pa)
    h0 = 0  # Base altitude (m)
    p = p0 * ((a + b * h) / (a + b * h0)) ** (-g / b / R)
    # Temperature
    h_ref = 0
    T_ref = T_sls
    dTdh = -6.5 / 1000  # K/m
    T = T_ref + dTdh * (h - h_ref)

    rho = p / R / T
    T = a + b * (h - h0)
    p_t = p * (1 + (gamma - 1) / 2 * M ** 2) ** (gamma / (gamma - 1))
    T_t = T * (1 + (gamma - 1) / 2 * M ** 2)

    return p, p_t, T, T_t, rho

# inputs:
# h: Altitude in meters
# V: airspeed
def isa_V(h, V):

    # International Standard Atmosphere
    # (dp/dh) = -(rho)*g
    # p = (rho)*R*T
    #
    # dp/p = -g*dh/(R*T)
    # T = a+b*h
    #
    # b != 0:
    # ln(p2/p1) = -g/(b*R)*ln((a+b*h2)/(a+b*h1))
    # p2 = p1*((a+b*h2)/(a+b*h1))^(-g/(b*R))
    #
    # b == 0:
    # ln(p2/p1) = -g*(h2-h1)/R/T
    # p2 = p1*exp(-g*(h2-h1)/R/T)
    
    g = 9.80665
    R = 287.058
    gamma = 1.4

    if h < 11000:
        # Troposphere
        a = 19 + 273 # Base temperature (K)
        b = -0.0065 # Temperature lapse rate (K/m)
        p0 = 108900 # Base pressure (Pa)
        h0 = -610 # Base altitude (m)
        p = p0 * ((a+b*h)/(a+b*h0))**(-g/b/R)
    elif h < 20000:
        # Tropopause
        a = -56.5 + 273 # Base temperature
        b = 0 # Temperature lapse rate
        p0 = 22632 # Base pressure
        h0 = 11000 # Base altitude
        p = p0 * exp(-g*(h-h0)/R/a)
    else:
        # Stratosphere (up to 32000 m)
        a = -56.5 + 273 # Base temperature
        b = 0.001 # Temperature lapse rate
        p0 = 5474.9 # Base pressure
        h0 = 20000 # Base altitude
        p = p0 * ((a+b*h)/(a+b*h0))**(-g/b/R)

    T = a + b * (h-h0)
        
    M = V/(gamma*R*T)**0.5
    p_t = p*(1+(gamma-1)/2*M**2)**(gamma/(gamma-1))
    T_t = T*(1+(gamma-1)/2*M**2)
    
    return p, p_t, T, T_t

if __name__ == "__main__":
    h_var = linspace(0, 12500, 200)  # [ft]
    fig, (ax1, ax3) = subplots(1, 2)
    # Reference 1959 ARDC Model
    h_ref = [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000]
    T_ref = [288.16, 281.66, 275.16, 268.67, 262.18, 255.69, 249.2, 242.71, 236.23, 229.74, 223.26, 216.78, 216.66, 216.66]
    p_ref = [101325, 89876, 79501, 70121, 61660, 54048, 47217, 41105, 35651, 30800, 26500, 22700, 19399, 16579]
    rho_ref = [1.225, 1.1117, 1.0066, 0.90926, 0.81935, 0.73643, 0.66011, 0.59002, 0.52578, 0.46706, 0.41351, 0.3648, 0.31194, 0.26659]
    T_ref_torr = [288.15, 281.65, 275.15, 268.65, 262.15, 255.65, 249.15, 242.65, 236.15, 229.65, 223.15, 216.65, 216.65,
             216.65]
    p_ref_torr = [101325, 89874, 79495, 70108, 61640, 54020, 47181, 41060, 35599, 30742, 26436, 22632, 19330, 16510]
    rho_ref_torr = [1.2250, 1.1117, 1.0065, 0.9091, 0.8191, 0.7361, 0.6597, 0.5895, 0.5252, 0.4663, 0.4127, 0.3639,
               0.3108, 0.2655]
    rho_plot = []
    p_plot = []
    T_plot = []
    p_t_plot = []
    T_t_plot = []
    for h in h_var:
        #p, p_t, T, T_t, rho = isa(h, 0.8)
        p, p_t, T, T_t, rho = troposphere(h, 0.8)
        rho_plot.append(rho)
        p_plot.append(p)
        T_plot.append(T)
        p_t_plot.append(p_t)
        T_t_plot.append(T_t)

    ax1.plot(h_var, T_plot, label='Troposphere')
    ax3.plot(h_var, rho_plot, label='Troposphere')

    rho_plot = []
    p_plot = []
    T_plot = []
    p_t_plot = []
    T_t_plot = []
    for h in h_var:
        # p, p_t, T, T_t, rho = isa(h, 0.8)
        p, p_t, T, T_t, rho = tropopause(h, 0.8)
        rho_plot.append(rho)
        p_plot.append(p)
        T_plot.append(T)
        p_t_plot.append(p_t)
        T_t_plot.append(T_t)

    ax1.plot(h_var, T_plot, label='Isotherm')
    ax1.plot(h_ref, T_ref, ls=':', label='1959 ARDC')
    ax1.plot(h_ref, T_ref_torr, ls=':', label='Torrenbeek')
    ax1.set_xlabel('Altitude [m]')
    ax1.set_ylabel('Temperature [K]')
    ax1.legend()
    ax3.plot(h_var, rho_plot, label='Tropopause')
    ax3.plot(h_ref, rho_ref, ls=':', label='1959 ARDC')
    ax3.plot(h_ref, rho_ref_torr, ls=':', label='Torrenbeek')
    ax3.set_xlabel('Altitude [m]')
    ax3.set_ylabel('Density [kg/m^3]')
    ax3.legend()
    show()

    h_var = linspace(0, 12500, 26)  # [ft]
    rho_plot = []
    p_plot = []
    T_plot = []
    p_t_plot = []
    T_t_plot = []
    for h in h_var:
        p, p_t, T, T_t, rho = isa(h, 0.8)
        rho_plot.append(rho)
        p_plot.append(p)
        T_plot.append(T)
        p_t_plot.append(p_t)
        T_t_plot.append(T_t)

    plot(h_var, T_plot, label='Calculated')
    plot(h_ref, T_ref, ls=':', label='1959 ARDC')
    plot(h_ref, T_ref_torr, ls=':', label='Torrenbeek')
    xlabel('Altitude [m]')
    ylabel('Temperature [K]')
    legend()
    show()

    plot(h_var, rho_plot, label='Calculated')
    plot(h_ref, rho_ref, ls=':', label='1959 ARDC')
    plot(h_ref, rho_ref_torr, ls=':', label='Torrenbeek')
    xlabel('Altitude [m]')
    ylabel('Density [kg/m^3]')
    legend()
    show()


    h_ref = [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000]
    rho_ref = [1.225, 1.1117, 1.0066, 0.90926, 0.81935, 0.73643, 0.66011, 0.59002, 0.52578, 0.46706, 0.41351, 0.3648, 0.31194, 0.26659]
    diff = []
    diff_torr = []
    for idx in range(0, len(h_ref)):
        p, p_t, T, T_t, rho = troposphere(h_ref[idx], 0.8)
        # p, p_t, T, T_t, rho = isa(h_ref[idx], 0.8)
        diff.append(rho / rho_ref[idx])
        diff_torr.append(rho / rho_ref_torr[idx])
    plot(h_ref, diff, label='1959 ARDC')
    plot(h_ref, diff_torr, label='Torrenbeek')
    xlabel('Altitude [m]')
    ylabel(r'Normalized Error $\rho / \rho_{ref}$')
    legend()
    show()


