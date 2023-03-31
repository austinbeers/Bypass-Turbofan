function  [Fsp, TSFC, eta_o, f, p7, T7, ue, p9, T9, uef] = turbofan(h, M, eta_d, eta_fn, eta_n, p_rf, BPR, eta_f, p_rc, eta_c, T04, r_b, eta_t, QR)
%         Calculate turbofan performance as a function of flight conditions and
%         cycle parameters.
%     
%         Inputs:
%             h       [m]         altitude
%             M                   flight Mach number
%             eta_d               inlet adiabatic efficiency
%             eta_fn              fan nozzle adiabatic efficiency
%             eta_n               core nozzle adiabatic efficiency
%             p_rf                fan stagnation pressure ratio p08/p02
%             BPR                 fan bypass ratio
%             eta_f               fan adiabatic efficiency
%             p_rc                overall compression ratio p03/p02
%             eta_c               compressor adiabatic efficiency
%             T04     [K]         burner exit stagnation temperature
%             r_b                 burner stagnation pressure ratio p04/p03
%             eta_t               turbine adiabatic efficiency
%             QR      [J/kg]      fuel heat release per unit mass
%     
%         Outputs:
%             Fsp     [N/(kg/s)]  specific thrust T / (\dot{m}_a)
%             TSFC    [lb/hr/lbf] thrust-specific fuel consumption

        R = 287.058;
        gamma_c = 1.4;
        gamma_h = 1.35;
    
        % free stream conditionsjump
        [pa, Ta, ~] = isa_func(h); 
        T0a = Ta * (1 + (gamma_c - 1) / 2 * M^2);
        u = M * (gamma_c * R * Ta)^0.5;
    
        % inlet/diffuser
        T02 = T0a;
        p02 = pa * (1 + eta_d * (T02/Ta - 1))^(gamma_c / (gamma_c - 1));
    
        % fan
        p08 = p02 * p_rf; 
        T08 = T02 * (1 + (p_rf^((gamma_c - 1) / gamma_c) - 1) / eta_f);
    
        % fan nozzle
        % First, assume unchoked
        T09 = T08;
        p9 = pa;
        p09 = p9 * (1 - eta_fn * (1 - (p9/p08)^((gamma_c - 1) / gamma_c)))^(-gamma_c / (gamma_c - 1));
        M9 = (2 / (gamma_c - 1) * ((p09/p9)^((gamma_c - 1) / gamma_c) - 1))^0.5;
        if M9 > 1  % nozzle must be choked, recalculate assuming M9 = 1
            fprintf('Fan nozzle is choked. \n')
            M9 = 1;
            p9 = p08 * (1 - (gamma_c - 1) / eta_fn / (gamma_c + 1))^(gamma_c / (gamma_c - 1));
            p09 = p9 / (1 + (gamma_c - 1) / 2)^(gamma_c / (gamma_c - 1));
        else
            fprintf('Fan nozzle is not choked. M_9 = %.2f \n', M9)
        end

        T9 = T09 / (1 + (gamma_c - 1) / 2 * M9^2);
        uef = M9 * (gamma_c * R * T9)^0.5;

        %  compressor
        p03 = p02 * p_rc;
        T03 = T02 * (1 + (p_rc^((gamma_c - 1) / gamma_c) - 1) / eta_c);
    
        % burner
        f = (T04 - T03) / (QR / (gamma_h * R / (gamma_h - 1)) - T04);
        p04 = p03 * r_b;
    
        % turbine
        T05 = T04 - (T03 - T02) - BPR * (T08 - T02);
        p05 = p04 * (1 - (1 - T05 / T04) / eta_t)^(gamma_h / (gamma_h - 1));
    
        % no afterburner
        T06 = T05;
        p06 = p05;
    
        % core nozzle
        % First, assume unchoked
        T07 = T06;
        p7 = pa;
        p07 = p7 * (1 - eta_n * (1 - (p7/p06)^((gamma_h - 1) / gamma_h)))^(-gamma_h / (gamma_h - 1));
        M7 = (2 / (gamma_h - 1) * ((p07/p7)^((gamma_h - 1) / gamma_h) - 1))^0.5;
        if M7 > 1  % nozzle must be choked, recalculate assuming M7 = 1
            fprintf('Core nozzle is choked. \n')
            M7 = 1;
            p7 = p06 * (1 - (gamma_h - 1) / eta_n / (gamma_h + 1))^(gamma_h / (gamma_h - 1));
            p07 = p7 / (1 + (gamma_h - 1) / 2)^(gamma_h / (gamma_h - 1));
        else
            fprintf('Core nozzle is not choked. M_7 = %.2f \n', M7)
        end

        T7 = T07 / (1 + (gamma_h - 1) / 2 * M7^2);
        ue = M7 * (gamma_h * R * T7)^0.5;
    
        % calculate performance metrics: specific thrust and TSFC
        Fsp = (1 + f) * ue + ue * (1 - pa / p7) / gamma_h / M7^2 + BPR * uef + BPR * uef * (1 - pa / p9) / gamma_c / M9^2 - (1 + BPR) * u;
        TSFC = f / Fsp;
        eta_o = u / QR / TSFC;
        TSFC = TSFC * 9.81 * 3600;  %  convert from kg/s/N to lb/hr/hp
    
        % output performance and nozzle exit states needed for area sizing
        end