% Fan Model 

% Constants
eta_pol0 = 0.90; 
m_tilde0 = 0.75;
a = 3.0;
delta_a = -0.5;
c = 3;
d = 6;
C = 2.5;
D = 15.0;

% Assume at design condition
p_tilde = 1.0; 
m_tilde = 1.0;

eta_pol = (1 - C * abs(p_tilde / (m_tilde^(a + delta_a - 1)) - m_tilde)^c - D * abs(m_tilde / m_tilde0 - 1)^d);

disp(eta_pol)
