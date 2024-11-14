clc
clear

%% Provided Skeleton Code
T = readtable(Filename);
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
theta2 = 30; % pitch
theta3 = 20; % theta3l

% DCM rotation based on a 3-2-1 rotation
C_BN = [cosd(theta2)*cosd(theta1), cosd(theta2)*sind(theta1), -sind(theta2);
        sind(theta3)*sind(theta2)*cosd(theta1)-cosd(theta3)*sind(theta1), sind(theta3)*sind(theta2)*sind(theta1)+cosd(theta3)*cosd(theta1), sind(theta3)*cosd(theta2);
        cosd(theta3)*sind(theta2)*cosd(theta1)+sind(theta3)*sind(theta1), cosd(theta3)*sind(theta2)*sind(theta1)-sind(theta3)*cosd(theta1), cosd(theta3)*cosd(theta2);
       ]

%% Functions
function isDCM = DCMcheck(A)
% This function checks if a matrix is a DCM
isDCM = true;

if (norm(A(1,:))~=1 || norm(A(2,:))~=1 || norm(A(3,:))~=1 ...
    || norm(A(:,1))~=1 || norm(A(:,1))~=1 || norm(A(:,1))~=1 )
    % Checks if any of the rows or columns are not equal to 1
    isDCM = false;

else if ((A*A') ~= eye(3))
    
end
% Checking if rows are unit vectors
A = A';
norm(A(1,:));
norm(A(2,:));
norm(A(3,:));


end