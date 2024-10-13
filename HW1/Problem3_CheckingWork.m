% AERSP 450 HW 1 Problem 3 (calculations only)
% Made by Nicholas Luis, PSU ID: 930841391

clear; clc;

% Constants
RE = 6378.137;
MU = 3.986 * (10^5);

% Given values
r1 = [14450.6, -1529.9, -6524];
r2 = [-6199.5, 14699.2, 8531.9];
p = 2.88*(RE);

h_dir = cross(r1,r2);
h_hat = h_dir./(norm(h_dir));

delta_f = acosd(sum(r1.*r2)/(norm(r1)*norm(r2)));

F = 1-norm(r2)/p * (1-cosd(delta_f));
G = norm(r2)*norm(r1)*sind(delta_f)/sqrt(MU*p);

v1 = (r2-F*r1)/G;

F_dot = (dot(r1,v1)/(p*norm(r1))) * (1-cosd(delta_f)) - (1/norm(r1))*sqrt(MU/p)*sin(delta_f);
G_dot = 1 - (norm(r1)/p) * (1-cosd(delta_f));

v2 = F_dot*r1 + G_dot*v1;

energy1 = 0.5*norm(v1)^2-MU/norm(r1); % Checks for constant energy
energy2 = 0.5*norm(v2)^2-MU/norm(r2); % Checks for constant energy
energy1 == energy2;

a = -0.5*MU*(0.5*norm(v1)^2-MU/norm(r1))^(-1);

e = sqrt(1-p/a);

h = cross(r1, v1);
e_vector = (1/MU)*(cross(v1, h)) - r1/(norm(r1));
e_hat = e_vector / e;

f1 = acosd(sum(e_hat.*(r1/norm(r1))));
f2 = acosd(sum(e_hat.*(r2/norm(r2))));