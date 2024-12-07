% AERSP 450 HW5 - Question III 
% Made by Nicholas Luis (PSU ID: 930841391)
% Remove semi-colons to see the output
clear; clc; close all;

%% Part A
% Known properties
h = 146.7; % mm
w = 71.5; % mm
t = 7.4; % mm
m = 189; % g

%Calculate moments of inertia 
Ix = (1/12) * m * (h^2 + t^2); 
Iy = (1/12) * m * (w^2 + t^2); 
Iz = (1/12) * m * (h^2 + w^2); 


%load('Experiment1.mat')
% load('Experiment2.mat')
load('Experiment3.mat')
wx = AngularVelocity.X;
wy = AngularVelocity.Y;
wz = AngularVelocity.Z;
T = AngularVelocity.Timestamp;
timeVec = datetime(T, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss.SSS');
t = seconds(timeVec - timeVec(1));

figure(1);
hold on;
    subplot(3,1,1);
    plot(t, wx);
%     title('Experiment 1 Data');
%     title('Experiment 2 Data');
    title('Experiment 3 Data');
    ylabel('\omega_x (rad/s)');
    %xline([5.102, 5.975], '--r'); % For experiment 1
    %     xline([3.18599, 4.11899], '--r'); % For experiment 2
    xline([4.6349, 5.42], '--r'); % For experiment 3

    subplot(3,1,2);
    plot(t, wy);
    ylabel('\omega_y (rad/s)');
%     xline([5.102, 5.975], '--r'); % For experiment 1
%     xline([3.18599, 4.11899], '--r'); % For experiment 2
    xline([4.6349, 5.42], '--r'); % For experiment 3

    subplot(3,1,3);
    plot(t, wz);
    xlabel('Time (s)');
    ylabel('\omega_z (rad/s)');
%     xline([5.102, 5.975], '--r'); % For experiment 1
%     xline([3.18599, 4.11899], '--r'); % For experiment 2
    xline([4.6349, 5.42], '--r'); % For experiment 3
hold off;

%% Part B
clear; clc; close all; 

% VALUES CONVERTED TO METERS AND KILOGRAMS TO GET JOULES DURING ENERGY CALCULATIONS
% Known properties
h = 146.7/1000; % m
w = 71.5/1000; % m
t = 7.4/1000; % m
m = 189/1000; % kg

%Calculate moments of inertia 
Ix = (1/12) * m * (h^2 + t^2); 
Iy = (1/12) * m * (w^2 + t^2); 
Iz = (1/12) * m * (h^2 + w^2); 

% Loading Data
load('Experiment1.mat')
wx1 = AngularVelocity.X;
wy1 = AngularVelocity.Y;
wz1 = AngularVelocity.Z;
T1 = AngularVelocity.Timestamp;
timeVec1 = datetime(T1, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss.SSS');
t1 = seconds(timeVec1 - timeVec1(1));

load('Experiment2.mat')
wx2 = AngularVelocity.X;
wy2 = AngularVelocity.Y;
wz2 = AngularVelocity.Z;
T2 = AngularVelocity.Timestamp;
timeVec2 = datetime(T2, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss.SSS');
t2 = seconds(timeVec2 - timeVec2(1));

load('Experiment3.mat')
wx3 = AngularVelocity.X;
wy3 = AngularVelocity.Y;
wz3 = AngularVelocity.Z;
T3 = AngularVelocity.Timestamp;
timeVec3 = datetime(T3, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss.SSS');
t3 = seconds(timeVec3 - timeVec3(1));

% Processing Data
KE1_hist = NaN(length(t1),1);
W1_hist = NaN(length(t1),1);
for i = 1:length(t1)
    KE1_hist(i) = (Ix*wx1(i)^2 + Iy*wy1(i)^2 + Iz*wz1(i)^2);
    W1_hist(i) = norm([Ix*wx1(i), Iy*wy1(i), Iz*wz1(i)]);
end

KE2_hist = NaN(length(t2),1);
W2_hist = NaN(length(t2),1);
for i = 1:length(t2)
    KE2_hist(i) = (Ix*wx2(i)^2 + Iy*wy2(i)^2 + Iz*wz2(i)^2);
    W2_hist(i) = norm([Ix*wx2(i), Iy*wy2(i), Iz*wz2(i)]);
end

KE3_hist = NaN(length(t3),1);
W3_hist = NaN(length(t3),1);
for i = 1:length(t3)
    KE3_hist(i) = (Ix*wx3(i)^2 + Iy*wy3(i)^2 + Iz*wz3(i)^2);
    W3_hist(i) = norm([Ix*wx3(i), Iy*wy3(i), Iz*wz3(i)]);
end

% Plotting 
figure(2)
hold on
plot(t1, KE1_hist, LineWidth=2)
plot(t2, KE2_hist, LineWidth=2)
plot(t3, KE3_hist, LineWidth=2)
legend('Experiment 1 ', 'Experiment 2', 'Experiment 3')
title('Total Energy')
xlabel("Time (s)")
ylabel("Energy (J)")
hold off
exportgraphics(gca,"HW5_Problem3_Energy.jpg");

figure(3)
hold on
plot(t1, W1_hist, LineWidth=2)
plot(t2, W2_hist, LineWidth=2)
plot(t3, W3_hist, LineWidth=2)
legend('Experiment 1 ', 'Experiment 2', 'Experiment 3')
title('Total Angular Momentum')
xlabel("Time (s)")
ylabel("Angular Momentum (kg*m^2/s)")
hold off
exportgraphics(gca,"HW5_Problem3_Momentum.jpg");





% The following are the regions of times that the phone is in the air (from part A)
times = [[5.102, 5.975];
         [3.18599, 4.11899];
         [4.6349, 5.42]];

% ---------Experiment 1 Progagation---------
t_range = times(1,:);

% Initial conditions
index = min(find(t1 >= 5.102));
wx1_0 = wx1(index);
wy1_0 = wy1(index);
wz1_0 = wz1(index);
w1_0 = [wx1_0, wy1_0, wz1_0]; % Vector of initial conditions for experiment 1

EOMs = @(t, omega) [ ((Iy - Iz) / Ix) * omega(2) * omega(3);
                     ((Iz - Ix) / Iy) * omega(3) * omega(1);
                     ((Ix - Iy) / Iz) * omega(1) * omega(2)];

% Progation using ode45
[t1_prop, omega1_prop] = ode45(EOMs, t_range, w1_0);

% Plotting
figure(4);
    subplot(3,1,1);
    plot(t1_prop, omega1_prop(:,1), '--r');
    hold on;
    title('Experiment 1 Propagated Data');
    plot(t1(index:index+length(t1_prop),:), wx1(index:index+length(t1_prop),:), 'b');
    hold off;
    legend('Analytical', 'Experimental')
    ylabel('\omega_x (rad/s)');

    subplot(3,1,2);
    plot(t1_prop, omega1_prop(:,2), '--r');
    hold on;
    plot(t1(index:index+length(t1_prop),:), wy1(index:index+length(t1_prop),:), 'b');
    hold off;
    legend('Analytical', 'Experimental')
    ylabel('\omega_x (rad/s)');

    subplot(3,1,3);
    plot(t1_prop, omega1_prop(:,3), '--r');
    hold on;
    plot(t1(index:index+length(t1_prop),:), wz1(index:index+length(t1_prop),:), 'b');
    hold off;
    legend('Analytical', 'Experimental')
    ylabel('\omega_x (rad/s)');
    xlabel('Time (s)');


    % ---------Experiment 2 Progagation---------
t_range = times(2,:);

% Initial conditions
index = min(find(t1 >= 3.18599));
wx2_0 = wx2(index);
wy2_0 = wy2(index);
wz2_0 = wz2(index);
w2_0 = [wx2_0, wy2_0, wz2_0]; % Vector of initial conditions for experiment 1

EOMs = @(t, omega) [ ((Iy - Iz) / Ix) * omega(2) * omega(3);
                     ((Iz - Ix) / Iy) * omega(3) * omega(1);
                     ((Ix - Iy) / Iz) * omega(1) * omega(2)];

% Progation using ode45
[t2_prop, omega2_prop] = ode45(EOMs, t_range, w2_0);

% Plotting
figure(5);
    subplot(3,1,1);
    plot(t2_prop, omega2_prop(:,1), '--r');
    hold on;
    title('Experiment 2 Propagated Data');
    plot(t2(index:index+length(t1_prop),:), wx2(index:index+length(t1_prop),:), 'b');
    hold off;
    legend('Analytical', 'Experimental')
    ylabel('\omega_x (rad/s)');

    subplot(3,1,2);
    plot(t2_prop, omega2_prop(:,2), '--r');
    hold on;
    plot(t2(index:index+length(t1_prop),:), wy2(index:index+length(t1_prop),:), 'b');
    hold off;
    legend('Analytical', 'Experimental')
    ylabel('\omega_x (rad/s)');

    subplot(3,1,3);
    plot(t2_prop, omega2_prop(:,3), '--r');
    hold on;
    plot(t2(index:index+length(t1_prop),:), wz2(index:index+length(t1_prop),:), 'b');
    hold off;
    legend('Analytical', 'Experimental')
    ylabel('\omega_x (rad/s)');
    xlabel('Time (s)');


        % ---------Experiment 3 Progagation---------
t_range = times(3,:);

% Initial conditions
index = min(find(t1 >= 4.6349));
wx3_0 = wx3(index);
wy3_0 = wy3(index);
wz3_0 = wz3(index);
w3_0 = [wx3_0, wy3_0, wz3_0]; % Vector of initial conditions for experiment 1

EOMs = @(t, omega) [ ((Iy - Iz) / Ix) * omega(2) * omega(3);
                     ((Iz - Ix) / Iy) * omega(3) * omega(1);
                     ((Ix - Iy) / Iz) * omega(1) * omega(2)];

% Progation using ode45
[t3_prop, omega3_prop] = ode45(EOMs, t_range, w3_0);

% Plotting
figure(6);
    subplot(3,1,1);
    plot(t3_prop, omega3_prop(:,1), '--r');
    hold on;
    title('Experiment 3 Propagated Data');
    plot(t3(index:index+length(t1_prop),:), wx3(index:index+length(t1_prop),:), 'b');
    hold off;
    legend('Analytical', 'Experimental')
    ylabel('\omega_x (rad/s)');

    subplot(3,1,2);
    plot(t3_prop, omega3_prop(:,3), '--r');
    hold on;
    plot(t3(index:index+length(t1_prop),:), wy3(index:index+length(t1_prop),:), 'b');
    hold off;
    legend('Analytical', 'Experimental')
    ylabel('\omega_x (rad/s)');

    subplot(3,1,3);
    plot(t3_prop, omega3_prop(:,2), '--r');
    hold on;
    plot(t3(index:index+length(t1_prop),:), wz3(index:index+length(t1_prop),:), 'b');
    hold off;
    legend('Analytical', 'Experimental')
    ylabel('\omega_x (rad/s)');
    xlabel('Time (s)');


%% Final Part

% Plotting ellipsoid
% Compute semi-axes of the ellipsoid from equation above
Figure(13);
hold on;
a =
b =
c =
% Create the ellipsoid
[x, y, z] = ellipsoid(0, 0, 0, a, b, c);
surf(x, y, z, 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor',...
'y', 'FaceAlpha', 0.5); % Uniform gray