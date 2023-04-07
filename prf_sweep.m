% Script to plot coefficient output 
prf_var = linspace(1.1, 1.8, 20);

alt = 10972.8;  % [m] 36,000 ft
mach = 0.78;  % 737 cruise mach
cfx_plot = zeros(1, length(prf_var));
cpx_plot = zeros(1, length(prf_var));

for i = 1:length(prf_var)
    [cpx, cfx] = calc_coefficients(alt, mach, prf_var(i));
    % Store values 
    cfx_plot(i) = cfx;
    cpx_plot(i) = cpx;
end

figure()
plot(cpx_plot, cfx_plot)
xlabel('Power Coefficient')
ylabel('Net Force Coefficient')