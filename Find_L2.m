% Code to find the L2 Lagrange point based on an initial (analytical) guess
% assuming small distance to the secondary.
% Author: Diogo Spínola
%
syms x
EarthMoon_mass = 6.0477e+24; % kg
Sun_mass = 1.9891e+30; %kg
mu = EarthMoon_mass/(Sun_mass + EarthMoon_mass)
L2e = 1 -mu + (mu./(3-2*mu))^(1/3)
f = @(x) x - (1 - mu)*(x + mu)/abs((x + mu)^3) - mu*(x - (1 - mu))/abs((x - (1 - mu))^3);
z = Newton(f,L2e,1e-15,1e15)