%% sensivity analysies for particle diameter in the perturbed CR3BP
d_range = 0.01:0.0005:0.1
L2pos = []
for d = d_range
    L2pos(end +1) = Find_L2(1,d);
end
%% Plotting
figure
plot(d_range*100, L2pos)
title('Sub-L_2 point as function of dust particle diameter.')
xlabel('Particle diameter [cm]')
ylabel('L_2 x adimensional coordinate')
saveas(gcf, 'Sub_L2_Sensivity.jpg')

%% Stability analysis - Unperturbed case
EarthMoon_mass = 6.0477e+24; % kg
Sun_mass = 1.9891e+30; %kg
mu = EarthMoon_mass/(Sun_mass + EarthMoon_mass);
x_L2 = Find_L2();
Uxx = -1 - 2*(1-mu)/(x_L2)^3 - 2*(mu)/(x_L2 - (1-mu))^3;
Uyy = -1 + (1-mu)/(x_L2)^3 + (mu)/(x_L2 - (1-mu))^3;
A = [
    0, 0, 1, 0
    0, 0, 0, 1
    -Uxx, 0, 0, 2
    0, -Uyy, -2, 0
    ];
e = eig(A)

x_L2_p = Find_L2(1);

d = 0.1; % Diameter of dust particle - [m]
rho = 1200; % Desity of the particle - [kg/m^3]
L_sun = 3.827*10^26; % Luminosity of the Sun - [J/s^1]
c = 299792458; % Speed of light in vacuum - [m/s]
R = 0.5; % Reflectivity of dust particle

mu_S = 1.327124421*10^20; % Sun's gravitational parameter - [m^3/s^2]

beta = (1+R)*(3*L_sun)/(8*pi*c*d*rho*mu_S); % Lightness number

Uxx_p = -1 - 2*(1-mu)*(1-beta)/(x_L2_p)^3 - 2*(mu)/(x_L2_p - (1-mu))^3;
Uyy_p = -1 + (1-mu)*(1-beta)/(x_L2_p)^3 + (mu)/(x_L2_p - (1-mu))^3;

A_p = [
    0, 0, 1, 0
    0, 0, 0, 1
    -Uxx_p, 0, 0, 2
    0, -Uyy_p, -2, 0
    ];
e_p = eig(A_p)