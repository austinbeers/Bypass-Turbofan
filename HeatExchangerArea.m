%% Heat Exchanger Model

M = 0.78; % Mach Number
rho = 0.38; % Density
T = 217; % Ambient temperature
gamma = 1.4; 
R = 287; % Gas Constant

a = sqrt(gamma*R*T); % Speed of sound
V = M*a; % Velocity

uh = 385; % Conductivity of copper
uc = 0.024; % Conductivity of air
uavg = (uh+uc)/2; % Average Conductivity 

cph = 2200; % Specific heat of glycol
cpc = 1004; % Specific heat of air

tc_in = T; % Inlet temperature
th_in = 394; % Temperature of working fluid

m_hot = 10; % Mass flow rate of working fluid
ch = m_hot*cph;

A = linspace(0.5,1.0,100); % Increasing area from 0.1 to 1.0
beta = linspace(75,1300,5);% heat transfer divided volume between plates

for j = 1:5
    for i = 1:100
        m_air(j,i) = rho*V*A(i); % Mass flow rate of incoming air

        cc(j,i) = m_air(j,i)*cpc;

        if cc(j,i) > ch
            cmax(j,i) = cc(j,i);
            cmin(j,i) = ch;
        else
            cmax(j,i) = ch;
            cmin(j,i) = cc(j,i);
        end

        cmin_cmax(j,i) = cmin(j,i)/cmax(j,i);

        NTU(j,i) = A(i)*uavg*beta(j)/cmin(j,i);

        e(j,i) = (1-exp(-NTU(j,i)*(1-(cmin_cmax(j,i)))))/(1-(cmin_cmax(j,i))*exp((-NTU(j,i))*(1-(cmin_cmax(j,i))))); % Efficiency
        tc_out(j,i) = ((e(j,i)*cmin(j,i)*(th_in - tc_in))/cc(j,i)) + tc_in; % Temperature output
        Q(j,i) = m_hot*(tc_out(j,i) - tc_in);
    end
    plot(A,Q(j,:))
    hold on;
end

title('Heating Addition vs. Inlet Area');
xlabel('Inlet Area');
ylabel('Heat Addition (Watts)');















