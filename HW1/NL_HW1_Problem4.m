% AERSP 450 HW 1 Problem 2
% Made by Nicholas Luis

%{
ISS (ZARYA)             
1 25544U 98067A   24259.04042691 -.00020782  00000+0 -36841-3 0  9993
2 25544  51.6359 230.2949 0007613 354.9391  85.5828 15.49088255472489

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