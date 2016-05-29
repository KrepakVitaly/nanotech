function [out, eta_fil, eta_ac, eta_hop, V_tunnel, current_full] = tafel( yi, V_app, apprx )
%TAFEL Summary of this function goes here
%   Detailed explanation goes here

e = 1.60217662e-19;
h = 6.62607004e-34;
kB = 1.38064852e-23;
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
m_eff = m_r * m_e;
C = 2.7;
T = 300 ;

fun = @(x) f(x, V_app, yi);
options = optimoptions('fsolve','Display','off');
x = fsolve(fun, apprx, options);
% disp(x);
eta_fil = x(1);
eta_ac = x(2);
eta_hop = x(3);
V_tunnel = x(4);

current_full = ((C * 3 * sqrt(2 * m_eff * ((V_tunnel/2)*e)) / 2 / yi * (e / h)^2 * ...
exp(- 4 * pi * yi / h * sqrt(2 * m_eff * ((V_tunnel/2)*e))) * A_fil*V_tunnel) ...
+ (j_0et*A_fil*(exp(-alpha*e*z*eta_fil/kB/T) - 1)));

current = j_0et * (exp((-alpha*e*z/kB/T*eta_fil)) - 1.0);
out = -M_me/z/e/rho_me * current;


end

