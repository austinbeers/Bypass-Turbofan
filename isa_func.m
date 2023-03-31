% Calculates thermodynamic state as a function of altitude for the International Standard Atmosphere.
% Inputs: h altitude in m Outputs: p pressure in Pa T temperature in K rho density in kg / m^3

function [p, T, rho] = isa_func(h)
    g = 9.80665; % m/s^2
    R = 287.058; % J/kg/K (dry air)
    p_sls = 101325; % Pa
    T_sls = 288.15; % K
    if h <= 11000
        h_ref = 0;
        p_ref = p_sls;
        T_ref = T_sls;
        dTdh = -6.5 / 1000; % K/m
        T = T_ref + dTdh * (h - h_ref);
        p = p_ref * (T / T_ref)^(-g / R / dTdh);
    elseif h <= 20000
        h_ref = 11000;
        [p_ref, T_ref, rho] = isa_func(h_ref);
        T = T_ref;
        p = p_ref * exp(-g * (h - h_ref) / R / T_ref);
    elseif h <= 32000
        h_ref = 20000;
        [p_ref, T_ref, rho] = isa_func(h_ref);
        dTdh = 1 / 1000;  % K/m
        T = T_ref + dTdh * (h - h_ref);
        p = p_ref * (T / T_ref)^(-g / R / dTdh);
    elseif h <= 47000
        h_ref = 32000;
        [p_ref, T_ref, ~] = isa_func(h_ref);
        dTdh = 2.8 / 1000; % K/m
        T = T_ref + dTdh * (h - h_ref);
        p = p_ref * (T / T_ref)^(-g / R / dTdh);
    end
    rho = p / R / T;
end
