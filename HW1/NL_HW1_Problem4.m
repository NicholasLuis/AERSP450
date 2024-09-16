% AERSP 450 HW 1 Problem 4
% Made by Nicholas Luis

%{

JUPITER 3 (ECHOSTAR 24) : This is the heaviest commericial satellite 
1 57479U 23108A   24257.28463005 -.00000150  00000+0  00000+0 0  9999
2 57479   0.0225 150.4751 0001204  19.0167 190.5107  1.00269368  4212

MOLNIYA 1-91            
1 25485U 98054A   24258.47584262 -.00000047  00000+0  00000+0 0  9999
2 25485  64.1500  21.4697 6821298 289.9022  12.2134  2.36441173198427

COSMOS 2569 (703K)      
1 57517U 23114A   24258.24908306  .00000071  00000+0  00000+0 0  9995
2 57517  64.9566  94.2948 0008739 295.2308  64.7601  2.13103887  8606

STARLETTE               
1 07646U 75010A   24258.93756172 -.00000102  00000+0  34095-4 0  9995
2 07646  49.8222 111.7558 0205657 107.9496 254.3932 13.82332614505978
%}

%{
To get position and velocity from orbital elements 
    Mean motion --> Period of orbit (T) --> (a)
    Eccentricity (e)
    (a) & (e) --> (r) & (v)
%}

clear; clc; close all;

% Delcaring constants
MU = 3.986 * (10^5); % km^3 / s^2
PI = 3.141592654;
%% ISS (ZARYA)             
% 1 25544U 98067A   24259.04042691 -.00020782  00000+0 -36841-3 0  9993
% 2 25544  51.6359 230.2949 0007613 354.9391  85.5828 15.49088255472489

% Orbit Data from the TLE
n = 15.49088255; % Revolutions (times 360 degrees) per day
e = 0.0007613; % Eccentricity
M0 = 85.5828; % Initial Mean Anomaly
w = 354.9391; % Argument of perigee
i = 51.6359; % Inclination
W = 230.2949; % Right ascension of ascending node
ERA = [0, 36, 10]; % GMST based on the provided website using epoch data

T = 86400/n; % Period of orbit (in seconds per 1 revolution)
a = ((T/(2*PI))^2 * MU)^(1/3);
t = linspace(0,T,1000); % Incremental time
t0 = t(1);
M = M0 + mod(360,n.*(t).*360); % Mean anomaly at every point in the orbit

% Converting M (mean anomaly) to E (eccentric anomaly) and f (true anomaly)
E = zeros(size(M));
for iter = 1 : length(M)
    E(iter) = mod(360,mean2true(M(iter), e));
end
E0 = E(1);
f0 = 2*atan( sqrt((1+e)/(1-e)) * tan(E(1)/2) ) * (180/PI);

% Getting initial radius and speed using derived equations
p = a*(1-e^2);
r0 = p / (1+e*cosd(f0));
v0 = sqrt(2*MU/r0 - MU/a);

% Calculating Langrange Coefficients using Block 2 (elliptic orbit)
F = 1 - (a/r0).*(1-cosd(E-E0)); 
G = (t-t0) - sqrt(a^3/MU).*((E-E0)*(PI/180)-sind(E-E0));

% Getting initial position and velocity vectors in perifocal frame
r = zeros(length(t), 3);
r(1,:) = [r0*cosd(f0), r0*sind(f0), 0]; % e, p, h hat directions
v = zeros(length(t), 3);
v(1,:) = [-sqrt(MU/p)*sind(f0), sqrt(MU/p)*(e+cosd(f0)), 0]; % e, p, h hat directions

% Calculating 1-orbit-worth of future positions
for iter = 2 : length(t) % Start at t = 2 because we already have initial conditions
    r(iter,:) = F(iter)*r(1,:) + G(iter)*v(1,:);
end

% Rotation Matrices
R3W = [cosd(W), -sind(W), 0;
       sind(W),  cosd(W), 0;
       0      ,  0      , 1 ]; % Rotation around axis 3 by W (Right ascension of ascending node) degrees
R1i = [1, 0      ,        0;
       0, cosd(i), -sind(i);
       0, sind(i),  cosd(i) ]; % Rotation around axis 1 by i (inclination) degrees
R3w = [cosd(w), -sind(w), 0;
       sind(w),  cosd(w), 0;
       0      ,  0      , 1 ]; % Rotation around axis 3 by w (argument of perigee) degrees

R_IP = R3w*(R1i*R3W);
ERA_degrees = [ERA(1)*15, ERA(2)*15/60, ERA(3)*15/60/60];
gamma0 = sum(ERA_degrees);

% Rotating position vectors to get it in ECEF frame (at every point in time)
for iter = 2 : length(t) % Start at t = 2 because we already have initial conditions
    g = gamma0 + (15/60)*t(iter); 
    R3g = [cosd(g), -sind(g), 0;
           sind(g),  cosd(g), 0;
           0      ,  0      , 1 ]; % Rotation around axis 3 by g (GMST angle) degrees
    R_EI = R_IP*R3g;
    r(iter,:) = R_EI*r(iter,:)';
end

% Converting position to lattitude and longitude
longitude = zeros(length(r));
latitude = zeros(length(r));
for iter = 1 : length(r)
    longitude(iter) = atan(r(iter, 2) / r(iter, 1)); % longtude equation using ECEF coords
    latitude(iter) = atan(r(iter, 3) / sqrt(r(iter, 1)^2 + r(iter, 2)^2 )); % lattitude equation using ECEF coords
end

%% Plots
% Eapen's Skeleton Code
figure
hold on
load('topo.mat','topo');
topoplot = [topo(:,181:360),topo(:,1:180)];
contour(-180:179,-90:89,topoplot,[0,0],'black');
grid on
grid minor
%%%%%%%%%
% WRITE CODE TO COMPUTE Longitude and latitude of the satellite in ECEF.
%%%%%%%%%
plot(longitude,latitude,'.','LineWidth',2)

%% Newton-Raphson Function
function E1 = mean2true(M, e) %Inputs a single number, Outputs a single number

    E0 = M; % Initial guess
    E1 = M + 1; % Ensures that loop exit condition is not prematurely called
    while (abs(E1-E0) > 0.00001) % Loops until the difference is negligible
        E0 = E1;
        f = E0 - e*sin(E0) - M;
        f_prime = 1 - e*sin(E0);
        E1 = E0 - (f / f_prime);
    end
end