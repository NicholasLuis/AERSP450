% This is the code for Aersp 450, HW 4, Question I
% Made by Nicholas Luis (PSU ID 930841391)

clc
clear

%% Provided Skeleton Code
T = readtable('SensorData.csv');
wx = T.wx;
wy = T.wy;
wz = T.wz;
% Step 1: Convert the time strings into datetime format
timeData = datetime(T.time, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''',...
'TimeZone', 'UTC');
% Step 2: Calculate time differences from the first time in the list
timeDifferences = timeData - timeData(1);
% Step 3: Convert the differences to seconds
t = seconds(timeDifferences);

%%

% Given euler angle rotations
theta1 = 30; % yaw
theta2 = 70; % pitch
theta3 = 20; % theta3l

% DCM rotation based on a 3-2-1 rotation
C_BN = [cosd(theta2)*cosd(theta1), cosd(theta2)*sind(theta1), -sind(theta2);
        sind(theta3)*sind(theta2)*cosd(theta1)-cosd(theta3)*sind(theta1), sind(theta3)*sind(theta2)*sind(theta1)+cosd(theta3)*cosd(theta1), sind(theta3)*cosd(theta2);
        cosd(theta3)*sind(theta2)*cosd(theta1)+sind(theta3)*sind(theta1), cosd(theta3)*sind(theta2)*sind(theta1)-sind(theta3)*cosd(theta1), cosd(theta3)*cosd(theta2);
       ]
DCMcheck(C_BN);

% Quaternion based on the DCM matrix using Sheppard Algo
Beta = SheppardAlgo(C_BN)

% Plotting the angular velocities as a function of time
figure(1)
hold on
plot(t, wx, LineWidth=2)
plot(t, wy, LineWidth=2)
plot(t, wz, LineWidth=2)
xlabel("Time (s)")
ylabel("Angular Velocity (deg / s)")
legend('wx', 'wy', 'wz')
hold off
exportgraphics(gca,"HW4_Problem1_AngVeloPlots.jpg");

%% Functions
function isDCM = DCMcheck(A)
% This function checks if a matrix is a DCM
isDCM = true;

% Checks if the rows and columns are unit vectors
for i = 1:3
    if (round(norm(A(i,:)), 10) ~= 1)
        isDCM = false;
        fprintf("The DCM is not valid!");
        return;
    end
    if (round(norm(A(:,i)), 10) ~= 1)
        isDCM = false;
        fprintf("The DCM is not valid!");
        return;
    end
end
% Checks if the DCM is orthonormal
if (round(A*A',10) ~= eye(3))
    isDCM = false;
    fprintf("The DCM is not valid!");
    return;
end
end


function BetaVec = SheppardAlgo(C)
% This funciton inputs some matrix C and does Sheppard's algorithm to
% compute the quaternion

% Equation 3.95
B0 = sqrt(0.25*(1+trace(C)));
B1 = sqrt(0.25*(1+2*C(1,1)-trace(C)));
B2 = sqrt(0.25*(1+2*C(2,2)-trace(C)));
B3 = sqrt(0.25*(1+2*C(3,3)-trace(C)));
BetaVec = [B0; B1; B2; B3];

% This part of Sheppard's algorithm leads to sum(B_i ^2) < 1 ???
%{
biggestB = max(BetaVec);

% Equation 3.96
if (B0 == biggestB)
    BetaVec(1) = 0.25*(C(2,3)-C(3,2))/B0;
    BetaVec(2) = 0.25*(C(3,1)-C(1,3))/B0;
    BetaVec(3) = 0.25*(C(1,2)-C(2,1))/B0;
    return;
elseif (B1 == biggestB)
    BetaVec(0) = 0.25*(C(2,3)-C(3,2))/B1;
    BetaVec(2) = 0.25*(C(1,2)+C(2,1))/B1;
    BetaVec(3) = 0.25*(C(3,1)+C(1,3))/B1;
    return;
elseif (B2 == biggestB)
    BetaVec(0) = 0.25*(C(3,1)-C(1,3))/B2;
    BetaVec(1) = 0.25*(C(1,2)+C(2,1))/B2;
    BetaVec(3) = 0.25*(C(2,3)+C(3,2))/B2;
    return;
elseif (B3 == biggestB)
    BetaVec(0) = 0.25*(C(1,2)-C(2,1))/B3;
    BetaVec(1) = 0.25*(C(3,1)+C(1,3))/B3;
    BetaVec(2) = 0.25*(C(2,3)+C(3,2))/B3;
    return;
end
%}
end
