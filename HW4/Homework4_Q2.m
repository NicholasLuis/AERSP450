% This is the code for Aersp 450, HW 4, Question II
% Made by Nicholas Luis (PSU ID 930841391)

clc
clear

%% Import Data
T = readtable('SensorData.csv');
wx = T.wx; % Roll rate
wy = T.wy; % Pitch rate
wz = T.wz; % Yaw rate
% Step 1: Convert the time strings into datetime format
timeData = datetime(T.time, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''',...
'TimeZone', 'UTC');
% Step 2: Calculate time differences from the first time in the list
timeDifferences = timeData - timeData(1);
% Step 3: Convert the differences to seconds
t = seconds(timeDifferences);

% Given euler angle rotations
theta1 = 30; % yaw
theta2 = 70; % pitch
theta3 = 20; % roll
% DCM rotation based on a 3-2-1 rotation
C_BN_original = [cosd(theta2)*cosd(theta1), cosd(theta2)*sind(theta1), -sind(theta2);
        sind(theta3)*sind(theta2)*cosd(theta1)-cosd(theta3)*sind(theta1), sind(theta3)*sind(theta2)*sind(theta1)+cosd(theta3)*cosd(theta1), sind(theta3)*cosd(theta2);
        cosd(theta3)*sind(theta2)*cosd(theta1)+sind(theta3)*sind(theta1), cosd(theta3)*sind(theta2)*sind(theta1)-sind(theta3)*cosd(theta1), cosd(theta3)*cosd(theta2);
       ];

%% Numerically Propogating the DCM at every timestep
C_BN = C_BN_original; % Creates a copy of the original DCM
C_hist_numerical = zeros(length(t)-1,9); % Stores the DCM at each timestep
C_BN_dot = zeros(3,3);

for i=1:length(t)-1
    % Saving the DCM as a single 9-element line
    C_hist_numerical(i,:) = mat2row(C_BN);
  
    % Propogating DCM
    omegaTilda = skewSymmetric([wx(i), wy(i), wz(i)]); % Creates the skew symmetric matrix of the angular velocity at a given time
    C_BN_dot = - omegaTilda * C_BN;
    C_BN = C_BN + (C_BN_dot*(t(i+1)-t(i)));
end

%% Analytically Propogating the DCM at every timestep
C_BN = C_BN_original; % Creates a copy of the original DCM
C_hist_analytic = zeros(length(t)-1,9); % Stores the DCM at each timestep

for j=1:length(t)-1
    C_hist_analytic(j,:) = mat2row(C_BN);
    omegaTilda = skewSymmetric([wx(j), wy(j), wz(j)]); % Creates the skew symmetric matrix of the angular velocity at a given time

    % Propagating DCM
    C_BN = expm( -omegaTilda.*(t(j+1)-t(j)) ) * C_BN;
end

%% Analytical vs Numerical Comparisons
error_hist = zeros(length(C_hist_numerical)-1,1); % Saves the error at every timestep
for i = 1:length(C_hist_numerical)

    % Gets the DCM at the given time step
    C_numeric = row2mat(C_hist_numerical(i,:));
    C_analyic = row2mat(C_hist_analytic(i,:));


    % Calculates error between two DCM's
    error_hist(i) = norm(C_numeric*C_analyic'-eye(3));
end

% Plotting
figure(1)
hold on
plot(t(1:length(t)-1), error_hist, LineWidth=2)
title('Time History of the Error' )
xlabel("Time (s)")
ylabel("Error magnitude")
hold off
exportgraphics(gca,"HW4_Problem2_ErrorPlot.jpg");

%% Yaw-Pitch-Roll vs Time
% Getting Euler Angle values from the DCMs
numericAngleHistory = zeros(length(t)-1, 3);
analyticAngleHistory = zeros(length(t)-1, 3);

for i = 1:length(t)-1
    % Numerical Data
    C_numeric = row2mat(C_hist_numerical(i,:));
    theta1 = atan2(C_numeric(2,1), C_numeric(1,1)); % Yaw
    theta2 = asin(C_numeric(3,1)); % Pitch
    theta3 = atan2(C_numeric(3,2), C_numeric(3,3)); % Roll
    numericAngleHistory(i,:) = [theta1, theta2, theta3]; % Saves yaw, pitch, and roll at this time step
    
    % Analytical Data
    C_analytic = row2mat(C_hist_analytic(i,:));
    theta1 = atan2(C_analytic(2,1), C_analytic(1,1)); % Yaw
    theta2 = asin(C_analytic(3,1)); % Pitch
    theta3 = atan2(C_analytic(3,2), C_analytic(3,3)); % Roll
    analyticAngleHistory(i,:) = [theta1, theta2, theta3]; % Saves yaw, pitch, and roll at this time step
end

% Plotting values
figure(2)
hold on
plot(t(1:length(t)-1), numericAngleHistory(:,1), LineWidth=2)
plot(t(1:length(t)-1), numericAngleHistory(:,2), LineWidth=2)
plot(t(1:length(t)-1), numericAngleHistory(:,3), LineWidth=2)
legend('Yaw ', 'Pitch', 'Roll')
title('Evolution of Yaw, Pitch, and Roll over Time (From Numerical DCM)' )
xlabel("Time (s)")
ylabel("Angle (Radians)")
hold off
exportgraphics(gca,"HW4_Problem2_EulerAnglesPlot.jpg");

figure(3)
hold on
plot(t(1:length(t)-1), analyticAngleHistory(:,1), LineWidth=2)
plot(t(1:length(t)-1), analyticAngleHistory(:,2), LineWidth=2)
plot(t(1:length(t)-1), analyticAngleHistory(:,3), LineWidth=2)
legend('Yaw ', 'Pitch', 'Roll')
title('Evolution of Yaw, Pitch, and Roll over Time (From Analytical DCM)' )
xlabel("Time (s)")
ylabel("Angle (Radians)")
hold off
exportgraphics(gca,"HW4_Problem2_EulerAnglesPlot2.jpg");
%% Functions
function matrixTilda = skewSymmetric(vec)
% This function inputs a vector and returns a skew symmetrix matrix
matrixTilda = [0, -vec(3), vec(2);...
               vec(3), 0, -vec(1);
               -vec(2), vec(1), 0;];
end

function rowMatrix = mat2row(inputMat)
% This function converts a 3x3 to a 9x1
rowMatrix = zeros(1,9);
ctr = 1;
    for iter = 1:3
        for jter = 1:3
            rowMatrix(ctr) = inputMat(iter,jter);
            ctr = ctr + 1;
        end
    end
end

function matrixOutput = row2mat(inputVec)
% This funciton inputs a 9x1 vector and converts it to a 3x3 matrix
matrixOutput = zeros(3,3);
ctr = 1;
    for jter = 1:3
        for iter = 1:3
            matrixOutput(iter,jter) = inputVec(ctr);
            ctr = ctr + 1;
        end
    end
end