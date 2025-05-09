function z = Find_L2(srp_active,d,rho)
arguments
    srp_active (1,1) double = 0;
    d (1,1) double = 0.1; % Diameter of dust particle - [m]
    rho (1,1) double = 1200; % Desity of the particle - [kg/m^3]
end
% Code to find the L2 Lagrange point based on an initial (analytical) guess
% assuming small distance to the secondary.
% Can also include spherical dust particle Solar Radiation Pressure
% Author: Diogo Spínola
%
syms x


% Solar radiation pressure related

L_sun = 3.827*10^26; % Luminosity of the Sun - [J/s^1]
c = 299792458; % Speed of light in vacuum - [m/s]
R = 0.5; % Reflectivity of dust particle

mu_S = 1.327124421*10^20; % Sun's gravitational parameter - [m^3/s^2]

beta = srp_active*(1+R)*(3*L_sun)/(8*pi*c*d*rho*mu_S); % Lightness number

% CR3BP related
EarthMoon_mass = 6.0477e+24; % kg
Sun_mass = 1.9891e+30; %kg
mu = EarthMoon_mass/(Sun_mass + EarthMoon_mass);

L2e = 1 -mu + (mu./(3-2*mu))^(1/3);
f = @(x) x - (1 - mu)*(x + mu)*(1-beta)/abs((x + mu)^3) - mu*(x - (1 - mu))/abs((x - (1 - mu))^3);
z = Newton(f,L2e,1e-15,1e15);
f(z)
end