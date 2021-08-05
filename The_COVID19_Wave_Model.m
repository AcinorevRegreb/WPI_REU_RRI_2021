% Final Code for The COVID-19 Wave Model (Berger, Tom, Yossefy)
% Model is based on the work of Gumel et al.
% Calls the function 'The_COVID19_Wave_Model_func' to solve the system

%% Data Sets
load('cumulative_deaths.mat')   % Data from Johns Hopkins
load('newcases.mat')            % Data from Johns Hopkins
load('hospitalizations.mat')    % Data from U.S. Department of Health & Human Services

% cumulative deaths from COVID-19 from 1/22/2020 to 7/22/2021
cumulativedeaths = [deaths20_0122to0229 deaths20_0301to0531 deaths20_0601to0930 deaths20_1001to0131 deaths21_0201to0531 deaths21_0601to0722];

% new cases per day of COVID-19 from 1/22/2020 to 7/22/2021
newcases = [newcases20_0122to0229 newcases20_0301to0531 newcases20_0601to0930 newcases20_1001to0131 newcases21_0201to0531 newcases21_0601to0722];

% current COVID-19 hospitalizations per day from 3/17/2020 to 3/7/2021
hospitalizations = [hosp20_0317to0531 hosp20_0601to0930 hosp20_1001to0131 hosp21_0201to0307];

%% Third Wave Data 

% for the purposes of our model, we will be fitting to third wave data
thirdwavecumulativedeaths = [deaths20_1001to0131 deaths21_0201to0531];
thirdwavenewcases = [newcases20_1001to0131 newcases21_0201to0531];

% to calculate daily deaths from cumulative death data, 
% we subtract day (i) cumulative deaths from day (i+1) cumulative deaths
thidwavedailydeaths = thirdwavecumulativedeaths(2:end)- thirdwavecumulativedeaths(1:end-1);

%% Initial Conditions

activecases=50000*14;       % # daily cases * average # of days infectious

S_u0=300e6;                 % Susceptible Unvaccinated (~US population)
S_v0=0;                     % Susceptible Vaccinated
E0=50000*2.5;               % Exposed (# daily cases * p.sigma)
P0=E0;                      % Presymptomatic
H0=30942;                   % Hospitalized (Department of Human & Health Services)
A0=0.3*(activecases-P0-H0); % Asymptomatic (30% of remaining active cases)
I0=0.7*(activecases-P0-H0); % Symptomatic (70% of remaining active cases)
R0=100000;                  % Recovered (assumed)
M0=thirdwavecumulativedeaths(1);    % Mortality (Johns Hopkins)

ICs = [S_u0;S_v0;E0;P0;A0;I0;H0;R0;M0];
%% Parameters for Third Wave Fit

% known parameters
p.pi = 1.2e4;       % natural birth rate (from Gumel paper)
p.mu = 1/(79*365);  % natural death rate (from Gumel paper)

p.xi_v = 0;         % vaccination rate per capita
p.epsilon_v = 0.95; % protective efficacy of vaccine

p.sigma = 1/2.5;    % incubation rate (E --> P)
p.alpha = 1/2.5;    % rate of (P --> I or A)
p.q = 0.7;          % proportion of pre-symp that become symptomatic
p.phi = 1/6;        % hospitalization rate for symptomatic

p.gamma_I = 1/10;   % recovery rate for infected
p.gamma_A = 1/5;    % recovery rate for asymptomatic
p.gamma_H = 1/8;    % recovery rate for hospitalized

% optimized parameters
p.beta_P0 = 5.35652305268018e-10;   % initial transmission rate of Presymptomatic class
p.beta_A0 = 1.32447536657924e-10;   % initial transmission rate of Asymptomatic class
p.beta_I0 = 9.80283810939473e-10;   % initial transmission rate of Symptomatic class
p.beta_H = 2.08273451259363e-22;    % fixed transmission rate of Hospitalized class
p.delta_I = 0.00172926752254379;    % disease-induced death rate of Symptomatic class
p.delta_H = 0.000325300035376653;   % disease-induced death rate of Hospitalized class

p.psi=0.00519468478446113;  % rate at which beta decreases after k reached
p.k=1352438.61962026;       % threshold in sigmoidal function
p.m=10;                     % width of sigmoidal function
p.p=0.02;                   % memory parameter

%% Solving the System of Equations

tspan=[0 242];     % number of days in between 10/1/2020 and 6/1/2021

% solving the system by using ode45 on The_COVID19_Wave_Model_func
[t,y]=ode45(@(t,x) (The_COVID19_Wave_Model_func(t,x,p)),tspan,[ICs;p.beta_I0]);

% solutions for populations
E = y(:,3);     % Exposed
P = y(:,4);     % Presymptomatic
A = y(:,5);     % Asymptomatic
I = y(:,6);     % Symptomatic
H = y(:,7);     % Hospitalized
M = y(:,9);     % Cumulative Mortality

% calculating daily deaths from the derivative of the Mortality equation (dM/dt)
dy=zeros(length(y),10);
for i=1:length(y)
     dy(i,:) = The_COVID19_Wave_Model_func(t(i), y(i,:), p);
end
dM=dy(:,9);
%% Plotting

% Cumulative Death Plot for 10/1/2020 to 6/1/2021
figure()
plot(0:242,thirdwavecumulativedeaths, 'o','markersize',4)   % Cumulative Death Data
hold on
plot(t,M,'linewidth',2)        % Cumulative Death from Simulation
title('The COVID-19 Wave Model: Cumulative Mortality Fit')
ylabel('Cumulative Disease-Induced Mortality') 
axis([0 242 2e5 6.5e5])
xticks([0 30 60 91 122 151 182 213 242])
xticklabels({'10/1/2020','11/1/2020','12/1/2020','1/1/2021','2/1/2021','3/1/2021','4/1/2021','5/1/2021','6/1/2021'})
xtickangle(45)
legend({'Data','COVID-19 Mortality from Simulation'},'Location','southeast')

% Daily Death Plot for 10/1/2020 to 6/1/2021
figure(2)
clf;
plot(1:242,thidwavedailydeaths,'o','markersize',4)  % Daily Death Data
hold on
plot(t,dM,'linewidth',2)       % Daily Death from Simulation
title('The COVID-19 Wave Model: Daily Mortality')
ylabel('Daily Disease-Induced Mortality') 
axis([0 242 0 5e3])
xticks([0 30 60 91 122 151 182 213 242])
xticklabels({'10/1/2020','11/1/2020','12/1/2020','1/1/2021','2/1/2021','3/1/2021','4/1/2021','5/1/2021','6/1/2021'})
xtickangle(45)
legend({'Data','COVID-19 Daily Mortality from Simulation'},'Location','northeast')

% New Case Plot for 10/1/2020 to 6/1/2021
figure(3)
clf;
plot(0:242,thirdwavenewcases,'o','markersize',4)   % New Case Data
hold on
plot(t,p.sigma*E,'linewidth',2)        % New Cases from Simulation 
title('The COVID-19 Wave Model: New Case Counts')
ylabel('Daily Case Counts of COVID-19') 
axis([0 242 1e4 6e5])
xticks([0 30 60 91 122 151 182 213 242])
xticklabels({'10/1/2020','11/1/2020','12/1/2020','1/1/2021','2/1/2021','3/1/2021','4/1/2021','5/1/2021','6/1/2021'})
xtickangle(45)
legend({'Data','Simulation of New Cases'},'Location','northeast')