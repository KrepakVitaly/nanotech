function F = f(x, V_app, x0)
%F Summary of this function goes here
%   Detailed explanation goes here

e = 1.60217662e-19;
h = 6.62607004e-34;
kB = 1.38064852e-23;
m_e = 9.10938356e-31;

% V_app = 0.15;

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
% x0 = L;

m_eff = m_r * m_e;

C = 2.7;
T = 300 ;

F(1) = kB * T / (1 - alpha) / e / z * log(A_fil/A_ac*(exp(- alpha * e * z / kB / T * x(1)) - 1) + 1) - x(2);

F(2) = x0*2*kB*T/a/z/e*asinh(j_0et/j_0hop*(exp(- alpha * e * z / kB / T * x(1)) - 1)) - x(3);

F(3) = x(2) - x(1) + x(3) - x(4);

F(4) = -V_app + ((C * 3 * sqrt(2 * m_eff * ((4+x(4)/2)*e)) / 2 / x0 * (e / h)^2 * ...
exp(- 4 * pi * x0 / h * sqrt(2 * m_eff * ((4+x(4)/2)*e))) * A_fil*x(4)) ...
+ (j_0et*A_fil*(exp(-alpha*e*z*x(1)/kB/T) - 1))) * (R_el + R_S + rho_fil*(L - x0) / A_fil) ...
+ x(4);
  
%         x(1) # eta_fil
%         x(2) # eta_ac
%         x(3) # eta_hop
%         x(4) # V_tunnel
end

