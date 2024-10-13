% AERSP 450 HW 1 Problem 2
% Made by Nicholas Luis, PSU ID: 930841391

clc; clear; close all;

% Constants
MU = 3.986 * (10^5); % km^3 / s^2
RE = 6378.137; % km, Radius of earth

% Known values
r = 2347 + RE; % km
v = 10.823; % km/s
f = 110;

%% Calculations
% Calulating orbital data
energy = (1/2)*v^2 - MU/r; 
a = -MU / (2*energy); % semimajor axis
h = r*v; % Angular momementum
e = sqrt(1 - h^2 / (MU*a)); % Eccentricity
H = 2*atanh(sqrt((e-1)/(e+1)) * tan(f/2)); % Hyperbolic anomaly
tf = (e * sinh(H) - H) * sqrt(abs(a^3) / MU); % Final time

% Calculating trajectory
ts = [0,tf];
s0 = [r;0;0; 0;v;0];
[t,s] = ode45(@(t,s) odefun1(t,s,mu), ts, s0);
xpos = flip(s(:,1));
ypos = flip(-s(:,2));
xpos = vertcat(xpos, s(:,1));
ypos = vertcat(ypos, s(:,2));

%% Plotting
% Plot the trajectory
fig1 = figure(1);
plt1 = plot(xpos,ypos,'r','LineWidth',3);
% Generate a circle (Earth)
Earth_iter = 0:1:360; % theta values from 0 to 360
xE = Re.*cosd(Earth_theta);
yE = Re.*sind(Earth_theta);
hold on;
pltE = plot(xE/1000,yE/1000, 'b','LineWidth',3);
grid on; axis equal;
legend([pltH,pltE], {'Hyperbolic Trajectory', 'Earth'});
title('Hyperbolic Trajectory of MESSENGER')
ylabel('y, [km]')
xlabel('x, [km]')
fontsize(fig, 'scale', 1.2)

%% Functions
% Note: It only iterates from 0-110. Because of symmetry, we can just
% mirror the plot
function dsdt = odefun1(t,s,mu)

 x = s(1);
 y = s(2);
 z = s(3);

 dxdt = s(4);
 dydt = s(5);
 dzdt = s(6);

 r = sqrt(x.^2 + y.^2 + z.^2);

 d2xdt2 = - mu .* x ./ r.^3;
 d2ydt2 = - mu .* y ./ r.^3;
 d2zdt2 = - mu .* z ./ r.^3;

 dsdt = NaN(6,1);
 dsdt(1) = dxdt;
 dsdt(2) = dydt;
 dsdt(3) = dzdt;
 dsdt(4) = d2xdt2;
 dsdt(5) = d2ydt2;
 dsdt(6) = d2zdt2;
end