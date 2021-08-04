% This contains the in host model main file

% load the data sets containing viral load data
% Found in Dogra et al., sourced from Vetter et al.
load('dataFromWholeBody.mat')
% Found in Esteban et al., sourced from Wolfel et al.
load('Fig5aDatapoints.mat')
load('Fig5bDataPoints.mat')

% scaling data from whole body
scaled_data = (viralload)/100;

% viral load values from fig 5a Esteban et al.
y_scaled = (transpose(10.^(y_average)));
% days for viral load data from fig 5a Esteban et al.
dpox = [4:22];

% viral load values from fig 5b Esteban et al.
y_final = (transpose(10.^y_final));
dpob = [3:20];
low_vl = (y_final)/100;
high_vl = (y_final)*100;

% viral load values from Vetter et al patient 2
day2 = [5,6,7,8,9,10,11,12];
viral_load2 = [2.5065e+10,3.2467e+09,1.7133e+09,9400000,32800000,3.5333e+08,7.1333e+05,20825];
viral_load2 = viral_load2./100;

% viral load values from Vetter et al patient 3
day3 = [6,7,8,9,10,11,12,13,24];
viral_load3 = [1.3475e+09,320000000,1.0600e+09,8.9388e+07,7.3333e+07,3.0667e+07,7.7333e+06,4200000,4.5620e+03];
viral_load3 = viral_load3./100;

%% Parameters
% Below are all the parameters used to create a fit for the Whole Body data
% mat, which is averaged data from Vetter et al

p.beta = 4.14;      % Infection rate of target cells
p.p_v = 18.72;      % Virus production rate in infected cells
p.eps = 1e-6;       % IFN effectiveness
p.alpha_pc = 4.04045058459382;    % influx of plasma from plasma to organ 
p.alpha_cp = 10.0886894165077;    % outflux of plasma from organ to plasma

p.d_t8 = 0.0139;    % Death rate of activated CD8 cells
p.d_C = 0.15;       % Cytopathic death rate of infected cells
p.d_I = 0.0037;     % Death rate of infected cells due to t8_hat (active)
p.d_Y = 0.0085;     % Elimination rate of viruses due to antibodies
p.d_F = 24;         % Degradation rate of IFN
p.d_A = 0.2;        % Death rate for activated APCs
p.d_t4 = 0.0185;    % Death rate of activated CD4 cells
p.d_S = 0.5;        % Death rate of short-lived plasma cells
p.d_L = 0.0083;     % Death rate of long-lived plasma cells

p.t4_gen = (10.^5.8); % Naïve CD4 steady state concentration
p.A_gen = (10.^5.3);  % Generic steady state concentration of APCs (diff vals per compartment)
p.B_gen = (1e5);      % Naïve B cells steady state concentration
p.t8_gen = (1e5);     % Naïve CD8 steady state concentration

p.g_F = 1.37;       % Baseline IFN production rate
p.p_F = 5.37;       % IFN production rate from infected cells
p.g_A = 0.1;        % Reequilibration rate of APCs
p.g_t8 = 0.0045;    % Reequilibration rate of naïve CD8 cells
p.g_B = 1.75;       % Reequilibration rate of naïve B cells
p.g_t4 = 0.0027;    % Reequilibration rate of naïve CD4 cells
p.g_S = 0.59;       % Production rate of antibodies from short-lived plasma cells
p.g_L = 1.34;       % Production rate of antibodies from long-lived plasma cells

p.t_B = 0.029*10^-3.5;     % Transition rate of naïve B cells into active B 
p.t_S = 564.8*10^-3.5;     % Transition rate of active B into short-lived plasma cells
p.t_L = 0.36*10^-3.5;      % Transition rate of active B into long-lived plasma cells
p.t_t = 0.0061*10^-3.5;    % Transition rate of naïve CD8 into active CD8 or CD4 into active CD4
p.t_A = 1.02e-9;           % Transition rate of naïve APCs into active APC 

p.lambda_Y = 0.0013; % Antibody loss rate due to viruses
p.c_Y = 0.04;        % Clearance rate of antibodies

p.E = 0.15;               % Hepatobiliary excretion rate
p.avg_lambda_M = 1.91;    % Viral elimination rate by macrophages in LRT and MPS

% these parameter values are for the low_vl data
%{
p.beta = 4;             % Infection rate of target cells
p.p_v = 12;             % Virus production rate in infected cells
p.eps = 1e-9;           % IFN effectiveness
p.t_B = 0.029*10^-2;    % Transition rate of naïve B cells into active B 
p.t_S = 564.8*10^-2;    % Transition rate of active B into short-lived plasma cells
p.t_L = 0.36*10^-2;     % Transition rate of active B into long-lived plasma cells
p.t_t = 0.0061*10^-2;   % Transition rate of naïve CD8 into active CD8 or CD4 into active CD4
p.t_A = 1e-9;           % Transition rate of naïve APCs into active APC 

% these parameter values are for the high_vl data

p.beta = 4;             % Infection rate of target cells
p.p_v = 168;            % Virus production rate in infected cells
p.eps = 2e-8;           % IFN effectiveness
p.t_B = 0.029*10^-3.7;  % Transition rate of naïve B cells into active B 
p.t_S = 564.8*10^-3.7;  % Transition rate of active B into short-lived plasma cells
p.t_L = 0.36*10^-3.7;   % Transition rate of active B into long-lived plasma cells
p.t_t = 0.0061*10^-3.7; % Transition rate of naïve CD8 into active CD8 or CD4 into active CD4
p.t_A = 9.7e-13;        % Transition rate of naïve APCs into active APC 

% these parameter values are for the individal patient data from Vetter et al., Patient 2

p.eps = 7.184e-08;      % IFN effectiveness
p.alpha_pc = 12.7;      % influx of plasma from plasma to organ 
p.alpha_cp = 0.9485;    % outflux of plasma from organ to plasma
p.t_B = 0.0054;         % Transition rate of naïve B cells into active B 
p.t_S = 1759.37;        % Transition rate of active B into short-lived plasma cells
p.t_L = 0.762;          % Transition rate of active B into long-lived plasma cells
p.t_t = 0.0025;         % Transition rate of naïve CD8 into active CD8 or CD4 into active CD4
p.t_A = 2.556e-14;      % Transition rate of naïve APCs into active APC
p.avg_lambda_M = 3;     % Viral elimination rate by macrophages in LRT and MPS

% these parameter values are for the individal patient data from Vetter et al., Patient 3

p.eps = 1e-08;          % IFN effectiveness
p.alpha_pc = 12.7;      % influx of plasma from plasma to organ 
p.alpha_cp = 0.9485;    % outflux of plasma from organ to plasma
p.t_B = 0.00494;        % Transition rate of naïve B cells into active B 
p.t_S = 2119.469;       % Transition rate of active B into short-lived plasma cells
p.t_L = 0.786;          % Transition rate of active B into long-lived plasma cells
p.t_t = 0.00204;        % Transition rate of naïve CD8 into active CD8 or CD4 into active CD4
p.t_A = 1.91e-14;       % Transition rate of naïve APCs into active APC
p.avg_lambda_M = 4.2;   % Viral elimination rate by macrophages in LRT and MPS
%}

%% Initial Conditions
% Initial conditions remain the same for all data sets and all fits
U0 = (10.^6.7)+(10.^6.7)+(10.^6.9)+(10.^6.6)+(10.^6.3)+(10.^6.9)+(10.^5.7);      % Uninfected
I0 = 100;                 % Infected
V_C0 = 1000;              % Viral load
F0 = 0;                   % IFN
V_P0 = 2500;              % Plasma
A0 = p.A_gen;             % Naive APC
A0_hat = 0;               % Active APC
T80 = p.t8_gen;           % Naive CD8
T80_hat = 0;              % Active CD8
B0 = p.B_gen;             % Naive B cell
B0_hat = 0;               % Active B cell
T40 = p.t4_gen;           % Naive CD4
T40_hat = 0;              % Active CD4
S_P0 = 0;                 % Short lived plasma
L_P0 = 0;                 % Long lived plasma
Y0 = 0;                   % Antibodies

y0 = [U0 I0 V_C0 F0 V_P0 A0 A0_hat T80 T80_hat B0 B0_hat T40 T40_hat S_P0 L_P0 Y0];

%% ODE Solver
% For the high viral load fit, the tspan becomes [0 50]
tspan = [0 26];
[t,y] = ode15s(@Our_In_Host_func2, tspan, y0, [], p);

%% Plotting

% Figure for viral load from Vetter et al. averaged
figure()
semilogy(t,y(:,3),'linewidth',2)      % V from simulation
hold on
semilogy(dpo,scaled_data, 'o','markersize',6) % viral load from Vetter et al.
legend({'Virus from Simulation', 'Viral Load Data (Vetter et al.)', ' Immune Response Data'})
xlabel('Time (days post infection)')
ylabel('Log10 Copies/mL')
title('Viral Load of SARS-CoV-2')
axis([0 26 1e0 1e9])

%{
% Figure for viral load from Vetter et al individual patients, Patient 2
figure()
semilogy(t,y(:,3),'linewidth',2)      % V from simulation
hold on
semilogy(day2,viral_load2, 'o','markersize',6)
legend({'Virus from Simulation', 'Viral Load Data'})
xlabel('Time (days post infection)')
ylabel('Log10 Copies/mL')
title('Patient 2 Viral Load of SARS-CoV-2')
axis([0 26 1e0 10^10])

% Figure for viral load from Vetter et al individual patients, Patient 3
figure()
semilogy(t,y(:,3),'linewidth',2)      % V from simulation
hold on
semilogy(day3,viral_load3, 'o','markersize',6)
legend({'Virus from Simulation', 'Viral Load Data', ' Immune Response Data'})
xlabel('Time (days post infection)')
ylabel('Log10 Copies/mL')
title('Patient 3 Viral Load of SARS-CoV-2')
axis([0 26 1e0 10^10])

% Figure for viral load from Wolfel et al. Patient A
figure()
semilogy(t,y(:,3),'linewidth',2)      % V from simulation
hold on
semilogy(dpox,y_scaled, 'o','markersize',6) % data from Wolfel et al. patient a
legend({'Virus from Simulation', 'Viral Load Data (Wolfel et al. Patient A)'})
xlabel('Time (days post infection)')
ylabel('Log10 Copies/mL')
title('Viral Load of SARS-CoV-2')
axis([0 26 1e0 1e9])

% Figure for simulating Low Viral Load from Wolfel et al.
figure()
semilogy(t,y(:,3),'linewidth',2)      % V from simulation
hold on
semilogy(dpob+5,low_vl, 'o','markersize',6) % data from Wolfel et al. patient b
legend({'Virus from Simulation', 'Viral Load Data (Low)'})
xlabel('Time (days post infection)')
ylabel('Log10 Copies/mL')
title('Viral Load of SARS-CoV-2')
axis([0 26 0 1e8])


% Figure for simulating High Viral Load from Wolfel et al.
figure()
semilogy(t,y(:,3),'linewidth',2)      % V from simulation
hold on
semilogy(dpob+5,high_vl, 'o','markersize',6) % data from Wolfel et al. patient b
legend({'Virus from Simulation', 'Viral Load Data (High)'})
xlabel('Time (days post infection)')
ylabel('Log10 Copies/mL')
title('Viral Load of SARS-CoV-2')
axis([0 50 0 1e10])


% Figure for Infected Cells over time
figure()
plot(t,y(:,2),'linewidth',2)
legend({'Infected Cells'})
title('Infected Cells')
xlabel('Time (days post infection)')
ylabel('Log10 Copies/mL')
axis([0 26 8e-17 3.3e7])

% Figure for Plasma over time
figure()
semilogy(t,y(:,5),'linewidth',2)  % plasma
legend({'Plasma'})
title('Plasma Compartment')
xlabel('Time (days post infection)')
ylabel('Log10 Copies/mL')
axis([0 26 -0.5e7 1e10])

% Figure for Naive CD4 over time
figure()
plot(t,y(:,12),'linewidth',2)
legend({'Naive CD4'})
title('Naive Cells')
xlabel('Time (days post infection)')
ylabel('Log10 Copies/mL')
axis([0 26 5.9e5 6.31e5])

% Figure for Naive APC over time 
figure()
plot(t,y(:,6),'linewidth',2)
legend({'Naive APC'})
title('Naive Cells')
xlabel('Time (days post infection)')
ylabel('Log10 Copies/mL')
axis([0 26 1.935e5 2e5])

% Figure for Naive CD8 and Naive B Cell over time
figure()
semilogy(t,y(:,8),'linewidth',2)
hold on
semilogy(t,y(:,10),'linewidth',2)
legend({'Naive CD8','Naive B Cell'})
title('Naive Cells')
xlabel('Time (days post infection)')
ylabel('Log10 Copies/mL')
axis([0 26 7.5e4 1.1e5])

% Figure for all Active Cells and Antibodies over time 
figure()
semilogy(t,y(:,7),'linewidth',2)  % active APCs
hold on
semilogy(t,y(:,9),'linewidth',2)  % active CD8s
semilogy(t,y(:,11),'linewidth',2) % active B Cells
semilogy(t,y(:,13),'linewidth',2) % active CD4
semilogy(t,y(:,16),'linewidth',2) % antibodies
axis([0 26 1e-8 1e10])
legend({'Active APC','Active CD8', 'Active Bcell','Active CD4','Antibodies'})
xlabel('Time (days post infection)')
ylabel('Log10 Copies/mL')
title('Active Cells')

% Figure for Short and Long Lived Plasma over time
figure()
hold on
semilogy(t,y(:,14),'linewidth',2) % short lived plasma
semilogy(t,y(:,15),'linewidth',2) % long lived plasma
legend({'Short-lived Plasma','Long-lived Plasma'})
title('Short and Long Lived Plasma')
xlabel('Time (days post infection)')
ylabel('Log10 Copies/mL')
axis([0 26 0 4580])

% Figure for IFN over time
figure()
semilogy(t,y(:,4),'linewidth',2)
legend({'IFN'})
title('IFN')
xlabel('Time (days post infection)')
ylabel('Log10 Copies/mL')
axis([0 26 1e-6 1e7])
%}