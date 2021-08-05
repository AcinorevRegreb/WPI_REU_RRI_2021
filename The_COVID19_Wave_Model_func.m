% Function for The COVID-19 Wave Model (Berger, Tom, Yossefy)
% Called by 'The_COVID19_Wave_Model'

function dydt = The_COVID19_Wave_Model_func(t,y,p)

% current sub-population sizes
S_u=y(1);
S_v=y(2);
E=y(3);
P=y(4);
A=y(5);
I=y(6);
H=y(7);
R=y(8);
M=y(9);

% current beta_I
beta_I=y(10);

dydt = zeros(10,1);

% updating presymp and asymp betas by normalizing beta_I and then multiplying by respective initial betas
% note: hospitalized beta is fixed, so does not need to be updated
beta_A=(beta_I/p.beta_I0)*p.beta_A0;
beta_P=(beta_I/p.beta_I0)*p.beta_P0;

% lambda: transmission from each infected class, update at beginning of time step
lambda = (beta_P * P + beta_I * I + beta_A * A + p.beta_H * H);    

% S_u: birth - transmission - vaxx - natural death
dydt(1) = p.pi + (-lambda - p.xi_v - p.mu) .* S_u; 

% S_v: vaxx - transmission - natural death
dydt(2) =  p.xi_v .* S_u - (1 - p.epsilon_v) * lambda .* S_v - p.mu .* S_v;          

% E: transmission - incubation - natural death
dydt(3) = lambda .* S_u + (1 - p.epsilon_v) * lambda .* S_v + (-p.sigma - p.mu) .* E;  

% P: incubation - symp/asymp - natural death
dydt(4) = p.sigma .* E + (-p.alpha - p.mu) .* P;    

% A: prop of pre-symp - recovery - natural death
dydt(5) = (1 - p.q)*p.alpha .* P + (-p.gamma_A - p.mu) .* A;

% I: prop of pre-symp - hospitalization - recovery - covid death - natural death
dydt(6) = p.q*p.alpha .* P + (-p.phi - p.gamma_I - p.delta_I - p.mu) .* I;  

% H: hospitalization - recovery - covid death - natural death
dydt(7) = p.phi .* I + (-p.gamma_H - p.delta_H - p.mu) .* H;

% R: symp recovery + asymp recovery + hosp recovery - natural death
dydt(8) = p.gamma_I .* I + p.gamma_A .* A + p.gamma_H .* H - p.mu .* R;

% M: symp covid death + hosp covid death
dydt(9) = p.delta_I .* I + p.delta_H .* H;

% beta_i: sigmoidal function + memory term
dydt(10) = -p.psi * ((H^p.m) / (p.k^p.m + H^p.m)) * beta_I + p.p * (p.beta_I0 - beta_I);