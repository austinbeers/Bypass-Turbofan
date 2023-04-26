function [mach] = calc_M(D, gamma)
    % Calculate mach number from corrected flow per unit area 
    mach_var = 0:0.01:1.0;
    D_m = mach_var ./ (1 + 0.5.*(gamma-1).*M.^2).^(0.5 * (gamma + 1)/(gamma - 1));
    
    % Interpolate mach from list of values 
    mach = interp1(D_m, mach_var, D);

end