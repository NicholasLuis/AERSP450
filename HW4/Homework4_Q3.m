% This is the code for Aersp 450, HW 4, Question II
% Made by Nicholas Luis (PSU ID 930841391)

clc
clear
close all

%% Import Data
T = readtable('SensorData.csv');
wx = T.wx; % Roll rate (omega 3)
wy = T.wy; % Pitch rate (omega 2)
wz = T.wz; % Yaw rate (omega 1)
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

% Initial Beta obtained (Copied from Part I)
Beta = [0.8050, -0.0088, 0.5824, 0.1126];

% Analytic Propagation (Copied from Part II)
C_BN = C_BN_original; % Creates a copy of the original DCM
C_hist_analytic = zeros(length(t)-1,9); % Stores the DCM at each timestep

for j=1:length(t)-1
    C_hist_analytic(j,:) = mat2row(C_BN);
    omegaTilda = skewSymmetric3x3([wx(j), wy(j), wz(j)]); % Creates the skew symmetric matrix of the angular velocity at a given time

    % Propagating DCM
    C_BN = expm( -omegaTilda.*(t(j+1)-t(j)) ) * C_BN;
end

% Getting Yaw, Pitch, and Roll values from the DCM (from Part II)
analyticAngleHistory = zeros(length(t)-1, 3);

for i = 1:length(t)-1    
    % Analytical Data
    C_analytic = row2mat(C_hist_analytic(i,:));
    theta1 = atan2(C_analytic(2,1), C_analytic(1,1)); % Yaw
    theta2 = asin(C_analytic(3,1)); % Pitch
    theta3 = atan2(C_analytic(3,2), C_analytic(3,3)); % Roll
    analyticAngleHistory(i,:) = [theta1, theta2, theta3]; % Saves yaw, pitch, and roll at this time step
end

%% Propagating the quaternion
B_hist = zeros(length(t)-1, 5); % Saves the quaternion vectors at every timestep (the last index is the magnitude)
B_hist(1,:) = [Beta, norm(Beta)];

for i = 1:length(t)-1
    B = skewSymmetric([wx(i),wy(i),wz(i)]);
    Phi = expm(0.5*B*(t(i+1)-t(i)));
    Bnew = Phi*B_hist(i,1:4)'; % Propagation
    B_hist(i+1, :) = [Bnew', norm(Bnew)];
end

% Checking if quaternion constraint is satisfied
constraintCheck = zeros(length(t), 1);
for j = 1:length(t)
    constraintCheck(j) = (B_hist(j,1)^2 + B_hist(j,2)^2 + B_hist(j,3)^2 + B_hist(j,4)^2);
end

% Plotting
figure(1)
hold on
plot(t(1:length(t)), round(constraintCheck,3), LineWidth=2)
title('Constraint check of Beta vectors' )
xlabel("Time (s)")
ylabel("Sum of (B_i^2)")
hold off
exportgraphics(gca,"HW4_Problem3_ConstraintCheck.jpg");

%% Converting Analytic Roll-Pitch-Yaw to Quaternion

% Converting Euler Angles to DCM at every timestep
DCM_hist = zeros(length(t)-1,9); % Creates a Nx9 matrix; Each 9-element row stores all the values of the DCM
for k = 1:length(t)-1
    DCM_hist(k,:) = mat2row( EulerAngles2DCM(analyticAngleHistory(k,:)) );
end

% Converting DCM to Quaternion at every timestep
B_Euler_hist = zeros(length(t)-1,4); % Each row in the matrix stores each of the 4 values of the quaternion
for L = 1:length(t)-1
    tempDCM = row2mat(DCM_hist(L,:)); % Convert the DCM row back into a 3x3
    B_Euler_hist(L,:) = DCM2Quaternion(tempDCM); % Converts the DCM into a Quaternion
end

%% Comparing the Analytic Quaternions to the Euler-Angle-Converted Quaternions
error = zeros(length(t)-1,1); % Vector to record the errors at each time step
for m = 1:length(t)-1
    deltaBeta = quaternionMultiplication(B_hist(m,1:4), B_Euler_hist(m,:));
    error(m) = norm(deltaBeta - [1,0,0,0]');
end

figure(2)
hold on
plot(t(1:length(t)-1), error, LineWidth=2)
title('Error History of the Quaternions obtained from Analyitcal vs Euler-Angle Methods')
xlabel("Time (s)")
ylabel("Magnitude of Error")
hold off
exportgraphics(gca,"HW4_Problem3_QuaternionErrorComparison.jpg");

%% Functions
function matrixTilda = skewSymmetric(vec)
% This function inputs a vector and returns a skew symmetrix matrix
matrixTilda = [0, -vec(1), -vec(2), -vec(3);...
               vec(1), 0, vec(3), -vec(2); ...
               vec(2), -vec(3), 0, vec(1); ...
               vec(3), vec(2), -vec(1), 0;];
end

function matrixTilda = skewSymmetric3x3(vec)
% This function inputs a vector and returns a skew symmetrix matrix
matrixTilda = [0, -vec(3), vec(2);...
               vec(3), 0, -vec(1);
               -vec(2), vec(1), 0;];
end

function DCM = EulerAngles2DCM(EulerAngles)
% This function converts Euler Angles to a DCM (from Part I.1)
theta1 = EulerAngles(1);
theta2 = EulerAngles(2);
theta3 = EulerAngles(3);
DCM = [cosd(theta2)*cosd(theta1), cosd(theta2)*sind(theta1), -sind(theta2);
        sind(theta3)*sind(theta2)*cosd(theta1)-cosd(theta3)*sind(theta1), sind(theta3)*sind(theta2)*sind(theta1)+cosd(theta3)*cosd(theta1), sind(theta3)*cosd(theta2);
        cosd(theta3)*sind(theta2)*cosd(theta1)+sind(theta3)*sind(theta1), cosd(theta3)*sind(theta2)*sind(theta1)-sind(theta3)*cosd(theta1), cosd(theta3)*cosd(theta2);
       ];
end

function BetaVec = DCM2Quaternion(C)
% This funciton inputs some 3x3 matrix C and computes the quaternion (from Part I.3)

% Equation 3.95
B0 = sqrt(0.25*(1+trace(C)));
B1 = sqrt(0.25*(1+2*C(1,1)-trace(C)));
B2 = sqrt(0.25*(1+2*C(2,2)-trace(C)));
B3 = sqrt(0.25*(1+2*C(3,3)-trace(C)));
BetaVec = [B0; B1; B2; B3];
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

function beta = quaternionMultiplication(B1, B2)
% This function performs quaternion multiplicaiton
B1_tilda = skewSymmetric(B1);
beta = B1_tilda*B2';
end