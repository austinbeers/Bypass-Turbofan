function [D] = calc_D(mach, gamma)
    % Calculate D(M) for corrected flow per unit area 
    D = mach / (1 + 0.5 * (gamma -1) * mach^2)^(0.5* ((gamma + 1) / (gamma - 1)));
end