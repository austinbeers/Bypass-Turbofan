from scipy.optimize import newton, root, root_scalar
from atmosphere import isa
from math import e
from utilities import calc_M, calc_D
from numpy import linspace, zeros, array, meshgrid, shape
from matplotlib.pyplot import plot, show, xlabel, ylabel, title, legend, gca, contour, contourf, colorbar

def hx(A, beta, m_dot_h, *args):
    # Extract args
    alt = args[0]  # m
    mach = args[1]


    # Constants
    R = 287  # [J/kg-K]
    cp = 1004  # [J/kg-K]
    gamma = 1.4

    # Ambient conditions
    pa, p0a, Ta, T0a, rho = isa(alt, mach)
    mdot = rho * mach * (gamma * R * Ta)**0.5 * A  # neglect choking

    # Heat exchanger design
    # Heat Addition
    Uh = 385  # [] Conductivity of copper
    Uc = 0.024  # [] Conductivity of air - should this be copper as well?
    Uavg = (Uh + Uc) / 2  # Average Conductivity

    Cph = 2200  # [J/kg-K] Specific heat of glycol
    Cpc = cp  # [J/kg-K] Specific heat of air

    Tc_in = Ta  # [K] Inlet temperature
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
    NTU = A * Uavg * beta / Cmin
    epsilon = (1 - e**(-1*NTU * (1 - Cmin / Cmax)))/(1 - (Cmin / Cmax) * e**(-1 * NTU * (1 - Cmin / Cmax)))
    q = epsilon * Ch * (Th_in - Tc_in)  # [W]

    return epsilon, q

if __name__ == "__main__":

    # Generate initial guess
    alt = 10972  # m
    mach = 0.78

    beta_var = linspace(100, 800, 40)  # [1/m^2]
    area_var = linspace(0.05, 0.3, 50)  # [m^2]

    area_plot, beta_plot = meshgrid(area_var, beta_var)
    epsilon_plot = zeros(shape(area_plot))
    q_plot = zeros(shape(area_plot))

    m_dot_h = 5  # [kg/s]
    for beta_idx in range(0, len(beta_var)):
        beta = beta_var[beta_idx]
        for area_idx in range(0, len(area_var)):
            area = area_var[area_idx]
            epsilon, q = hx(area, beta, m_dot_h, alt, mach)
            epsilon_plot[beta_idx][area_idx] = epsilon
            q_plot[beta_idx][area_idx] = q / 1000

    # Plot Contours
    contourf(area_plot, beta_plot, epsilon_plot)
    colorbar()
    xlabel(r'Area [$m^2$]')
    ylabel(r'$\beta$ [1/$m^2$]')
    title(r'Effectivenes, $\epsilon$, for $\dot{m}_h = 5 kg/s$')
    show()

    # Plot Contours
    contourf(area_plot, beta_plot, q_plot)
    colorbar()
    xlabel(r'Area [$m^2$]')
    ylabel(r'$\beta$ [1/$m^2$]')
    title(r'$\dot{Q}$ [kW] for $\dot{m}_h = 5 kg/s$')
    show()


    # Lower m_dot_h
    m_dot_h = 2  # [kg/s]
    for beta_idx in range(0, len(beta_var)):
        beta = beta_var[beta_idx]
        for area_idx in range(0, len(area_var)):
            area = area_var[area_idx]
            epsilon, q = hx(area, beta, m_dot_h, alt, mach)
            epsilon_plot[beta_idx][area_idx] = epsilon
            q_plot[beta_idx][area_idx] = q / 1000

    # Plot Contours
    contourf(area_plot, beta_plot, epsilon_plot)
    colorbar()
    xlabel(r'Area [$m^2$]')
    ylabel(r'$\beta$ [1/$m^2$]')
    title(r'Effectivenes, $\epsilon$, for $\dot{m}_h = 2 kg/s$')
    show()

    # Plot Contours
    contourf(area_plot, beta_plot, q_plot)
    colorbar()
    xlabel(r'Area [$m^2$]')
    ylabel(r'$\beta$ [1/$m^2$]')
    title(r'$\dot{Q}$ [kW] for $\dot{m}_h = 2 kg/s$')
    show()