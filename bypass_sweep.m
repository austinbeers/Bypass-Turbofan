% run turbofan script 
% This part of the code calls the turbofan() function and ... when the script is called from the terminal or Python shell wit
% Note: numpy and matplotlib Python packages must be installed, otherwise the script won't run

hm = 10668;
fprintf('Cruise \n')
[Fsp, TSFC, eta_o, f, p7, T7, u7, p9, T9, u9] = turbofan(hm, 0.8, 0.99, 0.99, 0.99, 1.51, 14, 0.94, 45, 0.91, 1700, 1.00, 0.95, 43e6);

fprintf('eta_o = %.2f \n', eta_o)
fprintf('TSFC = %.3f lb/hr/lbf \n', TSFC)
fprintf('F/\dot{m}_a = %.2f N/(kg/s) \n', Fsp)
fprintf('Nozzle sizing \n')
T = 32928; % 7400 lbf in Newtons
mdot_a = T / Fsp; % get core mass flow from thrust and Fsp
fprintf('Core mass flow: %.2f kg/s \n', mdot_a)
% core nozzle mass flow, density, velocity --> area
mdot_7 = mdot_a * (1 + f);
rho7 = p7 / 287.058 / T7;
A7 = mdot_7 / rho7 / u7;
fprintf('A_7 = %.1f sq m \n', A7)
% bypass nozzle mass flow, density, velocity --> area
mdot_9 = mdot_a * 14;
rho9 = p9 / 287.058 / T9;
A9 = mdot_9 / rho9 / u9;
fprintf('A_9 = %.1f sq m \n', A9)
% Part 3
fprintf('Sea Level Static \n')
% Loop over range of bypass ratios to see where specified = calculated
BPR_spec = linspace(0, 10.1, 200);
BPR_calc = zeros(1, length(BPR_spec));
for ii = 1:length(BPR_spec)
    [~, ~, ~, f, p7, T7, u7, p9, T9, u9] = turbofan(0, 0, 0.99, 0.99, 0.99, 1.51, BPR_spec(ii), 0.94, 45, 0.91, 1700, 1.00, 0.95, 43e6);
    rho_7 = p7 / 287.058 / T7;
    rho_9 = p9 / 287.058 / T9;
    mdot_7 = rho_7 * u7 * A7;
    mdot_9 = rho_9 * u9 * A7;
    BPR_calc(ii) = mdot_9 / (mdot_7 / (1 + f));
end
% plotting
figure()
plot(BPR_spec, BPR_spec); hold all; 
plot(BPR_spec, BPR_calc)
plot(0.419, 0.419)
xlabel('Specified bypass ratio')
ylabel('Calculated bypass ratio')
axis([0, 10, 0, 10])
axis equal
%  calcluate force for BPR = 0.419 design
[Fsp, TSFC, eta_o, f, p7, T7, u7, p9, T9, u9] = turbofan(0, 0, 0.99, 0.99, 0.99, 1.51, 0.419, 0.94, 45, 0.91, 1700, 1.00, 0.95, 43e6);
mdota = p7 * u7 * A7 / 287.058 / T7;
F = Fsp * mdota / 9.81 * 2.2;
fprintf('T_SLS = %.0f lbf\n', F)