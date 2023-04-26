% Calculate streamwise force and power coefficients
function [cp, cfx, station_mach] = calc_coefficients(alt, M, prf, Ad, An)
    % Air
    gamma = 1.4;
    R = 287.05;  % J/ kg-K

    % free stream conditions
    [pa, Ta, rho] = isa_func(alt);
    T0a = Ta * (1 + (gamma - 1) / 2 * M^2);
    p0a = pa * (1 + (gamma - 1) / 2 * M^2)^(gamma / (gamma - 1));
    u = M * (gamma * R * Ta)^0.5;
    
    % isentropic diffuser 
    T01 = T0a;
    p01 = p0a;
    
    % fan pressure ratio
    p02 = p01 * prf;
    T02 = T01 * (1 + (prf^((gamma - 1) / gamma) - 1));
    
    % isentropic nozzle
    T0e = T02;
    p0e = p02;
    pe = pa;
    Me = ((2 / (gamma - 1)) * ((p0e / pe)^((gamma - 1) / gamma) - 1))^0.5;
    if Me > 1  % check if nozzle is choked
        fprintf('Fan nozzle is choked. \n')
        Me = 1;
        pe = p0e * (1 + ((gamma - 1) / 2) * Me^2)^(-1 * gamma / (gamma -1));
    else
        fprintf('Fan nozzle is not choked. M_e = %.2f \n', Me)
    end
    Te = T0e / (1 + ((gamma - 1) / 2) * Me^2);
    ue = Me * (gamma * R  * Te)^0.5;
    De = Me / ((1 + ((gamma - 1) / 2) * Me^2)^(0.5 * (gamma + 1) / (gamma - 1)));

    cp = ((De * p0e * gamma^0.5) / (R * T0e)^0.5) * (gamma * R / (gamma - 1)) * T01 * (prf - 1)^((gamma - 1) / gamma) / (0.5 * rho * u^3);
    cfx = (((De * p0e * gamma^0.5) / (R * T0e)^0.5) * (ue - u) + (pe - pa)) / (0.5 * rho * u^2);
    
    % Calculate Mach and Areas at each station
    % Station 2
    D2 = De * An;
    M2 = calc_M(D2, gamma);

    % Station 1 
    M1 = M2;
    D1 = calc_D(M1, gamma);

    % Station Inlet
    Di = D1 * Ad;
    Mi = calc_M(Di, gamma);

    station_mach = [Mi, M1, M2, Me];
end

