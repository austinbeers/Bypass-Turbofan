% Density Compressible Flow Relation

gamma = 1.4;
mach_var = linspace(0, 1.0, 50);
rho = zeros(1, length(mach_var));

for i = 1:length(mach_var)
    mach = mach_var(i);
    rho(i) = 1 / ((1 + 0.5 * (gamma - 1) * mach^2)^(1 / (gamma - 1)));
end

figure()
plot(mach_var, rho)
xlabel('Mach')
ylabel('rho / rho_t')