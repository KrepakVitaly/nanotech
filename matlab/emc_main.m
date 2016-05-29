clear;
disp('PARTY IS ON!!!!')
e = 1.60217662e-19;
h = 6.62607004e-34;
kB = 1.38064852e-23;
m_e = 9.10938356e-31;

M_me       = 1.79e-25  ; % kg
z          = 1.0  ; % n/a
rho_me     = 10.49e3  ; % kg m-3
m_r        = 0.023  ; % n/a
alpha      = 0.3  ; % n/a
j_0et      = 3.2e5  ; % A m-2
DeltaG_et  = 0.6 * e  ; % Joule
j_0hop     = 1.1e11  ; % A m-2
a          = 0.25e-9  ; % m
DeltaG_hop = 0.32 * e  ; % Joule
DeltaG_nuc = 0.80 * e  ; % Joule
t_0nuc     = 2e-8  ; % s
N_c        = 3.0  ; % n/a
A_ac       = 804.25e-18  ; % m2
A_fil      = 12.57e-18  ; % m2
A_is       = 12.57e-18  ; % m2
L          = 30.0e-9  ; % m
rho_fil    = 1.7e-8  ; % Ohm m
R_el       = 76.4  ; % Ohm
R_S        = 1.0e6  ; % Ohm
x0 = L;
m_eff = m_r * m_e;
C = 2.7;
T = 300 ;

V_app = 2;

% fun = @(x) f(x, V_app, x0);
% x0 = [0,0,0,0];
% options = optimoptions('fsolve','Display','off');
% [x,fval,exitflag,output] = fsolve(fun,x0, options);
% disp(x);

% eta=x(4);

% t_nuc = t_0nuc * exp(DeltaG_nuc/kB/T) * exp(-(N_c + alpha)*z * e * eta/kB/T);
        

% x = ones(1,10);
time = 5e-8;
n_steps = 500;
step = time/n_steps;

y = zeros(n_steps+1, 1);
y_s = zeros(n_steps+1, 1);
eta_fil = zeros(n_steps+1, 1);
eta_ac = zeros(n_steps+1, 1);
eta_hop = zeros(n_steps+1, 1);
V_tunnel = zeros(n_steps+1, 1);
current_full = zeros(n_steps+1, 1);

y(1) = L;
yi = 0.0;
apprx = [0,0,0,0];
for i = 1:n_steps
    disp(i)
    [yi, eta_fil(i), eta_ac(i), eta_hop(i), V_tunnel(i), current_full(i)] = tafel(y(i), V_app, apprx);
    y_s(i+1) = y(i) + step * yi;
    [ysi1, eta_fil(i), eta_ac(i), eta_hop(i), V_tunnel(i), current_full(i)] = tafel(y_s(i+1), V_app, apprx);
    y(i+1) = y(i) + step/2 * (yi + ysi1);
end

plot(eta_ac);

