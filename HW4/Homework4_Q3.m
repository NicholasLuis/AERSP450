% This is the code for Aersp 450, HW 4, Question II
% Made by Nicholas Luis (PSU ID 930841391)

clc
clear

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
C_BN = [cosd(theta2)*cosd(theta1), cosd(theta2)*sind(theta1), -sind(theta2);
        sind(theta3)*sind(theta2)*cosd(theta1)-cosd(theta3)*sind(theta1), sind(theta3)*sind(theta2)*sind(theta1)+cosd(theta3)*cosd(theta1), sind(theta3)*cosd(theta2);
        cosd(theta3)*sind(theta2)*cosd(theta1)+sind(theta3)*sind(theta1), cosd(theta3)*sind(theta2)*sind(theta1)-sind(theta3)*cosd(theta1), cosd(theta3)*cosd(theta2);
       ];

%% Propagating the quaternion


%% Functions
function matrixTilda = skewSymmetric(vec)
% This function inputs a vector and returns a skew symmetrix matrix
matrixTilda = [0, -vec(1), -vec(2), -vec(3);...
               vec(1), 0, vec(3), -vec(2); ...
               vec(2), -vec(3), 0, vec(1); ...
               vec(3), vec(2), -vec(1), 0;];
end