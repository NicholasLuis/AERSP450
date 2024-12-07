% AERSP 450 HW5 - Question 1 
% Made by Nicholas Luis (PSU ID: 930841391)
% Remove semi-colons to see the output

%% Parts A & B - Basic Data
clc; clear; close all;

% Step 1
grayImg = imread('StarField.png');
imshow(grayImg);

% Step 2
binaryImg = imbinarize(grayImg, 'adaptive'); %Adaptive thresholding
imshow(binaryImg)

% Step 3
stats = regionprops(binaryImg, 'Centroid');
centroids = cat(1, stats.Centroid); %Extract x and y pixel coordinates

% Step 4 : Comparing extracted coords to manually-obtained coords
Star1 = [86,19];
Star2 = [150,14];
Star3 = [10,114];
Star4 = [60,93];
Star5 = [150,101];
starCoords = [Star1; Star2; Star3; Star4; Star5];

Cents = NaN(length(starCoords),2);
for i=1:size(starCoords)
    [M,I] = min(((centroids(:,1)-starCoords(i,1)).^2 + (centroids(:,2)-starCoords(i,2)).^2).^(0.5));
    Cents(i,:) = centroids(I,:); % Saves new centroid coordinates
end
clear i;

% Step 5: Verifying coordinates
imshow("StarField.png")
hold on
plot(centroids(:,1),centroids(:,2),'.r')
hold on
plot(Cents(:,1),Cents(:,2),'og')

% Step 6: Converting pixel coords into spatial
[L,W] = size(binaryImg); % Gets image dimensions
P = (Cents - W/2) * tan((4*pi/180)/2);
P = cat(2,P,ones(5,1));
for j = 1:length(P)
    P(j,:) = P(j,:)./norm(P(j,:)); % Normalize vectors to get unit vector
end
clear j;

% (Part B)
StarI1 = [-0.921069884293268 -0.342599924017704 0.185082577005643];
StarI2= [-0.476980639452282 -0.711962047538100 0.515363476056510];
StarI3= [-0.496592767571065 0.816782636490539 -0.293703503424244];
StarI4= [-0.736236774114155 0.676644212846351 0.010393347079969];
StarI5 = [ 0.304668730748568 -0.299446239108458 0.904161995655567];
starCoordsInertial = [StarI1; StarI2; StarI3; StarI4; StarI5];

close all;
%% Part C (TRIAD Method)
% Step 1: Computing DCM
b1 = P(2,:)';
b2 = P(3,:)';
b3 = cross(b1,b2) / norm(cross(b1,b2)); 

r1 = starCoordsInertial(2,:)';
r2 = starCoordsInertial(3,:)';
r3 = cross(r1,r2) / norm(cross(r1,r2));

B = [b1, b3, cross(b1,b3)];
R = [r1, r3, cross(r1,r3)];

C_BN_triad = B*R';

% Step 2
errorMat_Triad = NaN(3,length(starCoordsInertial));
for k = 1:length(starCoordsInertial)
    errorMat_Triad(:,k) = P(k,:)' - (C_BN_triad*starCoordsInertial(k,:)');
end
clear k;
errorTriad = norm([errorMat_Triad(1,:), errorMat_Triad(2,:), errorMat_Triad(3,:)]);

%% Part D (OLAE Method)
S_matrix = NaN(3*length(P), 3);
D_vector = NaN(length(P), 3);

for n = 1:length(P)
    S = P(n,:)' + starCoordsInertial(n,:)';
    D_vector(n,:) = P(n,:)' - starCoordsInertial(n,:)';

    cntr = 1 + 3*(n-1);
    S_matrix(cntr:3*n, :) = skewSymmetric(S);
end
clear n; clear cntr;

D_long = reshape(D_vector',[],1); % Long vector

q = pinv(S_matrix)*D_long;

% DCM relating body to inertial frame using OLAE method
C_BN_olae = ((eye(3) + skewSymmetric(q)) / (eye(3) - skewSymmetric(q)))';

% Calculating error
errorMat_OLAE = NaN(3,length(starCoordsInertial));
for k = 1:length(starCoordsInertial)
    errorMat_OLAE(:,k) = P(k,:)' - (C_BN_olae*starCoordsInertial(k,:)');
end
clear k;
errorOlae = norm([errorMat_OLAE(1,:), errorMat_OLAE(2,:), errorMat_OLAE(3,:)]);


%% Part E (Davenport q-method)
B_mat = zeros(3,3);
for i = 1:length(P)
    B_mat = B_mat + ( P(i,:)' * starCoordsInertial(i, :) );
end

sigma = trace(B_mat);
S = B_mat + B_mat';
Z = [B_mat(2,3) - B_mat(3,2);
     B_mat(3,1) - B_mat(1,3);
     B_mat(1,2) - B_mat(2,1); ];

% K matrix
K = [sigma, Z';
     Z(1), S(1,:) - [sigma, 0, 0];
     Z(2), S(2,:) - [0, sigma, 0];
     Z(3), S(3,:) - [0, 0, sigma]; ];

% Eigen-stuff
[eVec, eVal] = eig(K);
[~, I] = max(diag(eVal));
q = eVec(:,I);

% Converting quaternion to DCM
C_BN_dave = [ q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2, 2*(q(2)*q(3) + q(1)*q(4)), 2*(q(2)*q(4) - q(1)*q(3)) ;
              2*(q(2)*q(3) - q(1)*q(4)), q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2, 2*(q(3)*q(4) + q(1)*q(2)) ;
              2*(q(2)*q(4) + q(1)*q(3)), 2*(q(3)*q(4) - q(1)*q(2)), q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2 ;]

% Calculating error
errorMat_dave = NaN(3,length(starCoordsInertial));
for k = 1:length(starCoordsInertial)
    errorMat_dave(:,k) = P(k,:)' - (C_BN_dave*starCoordsInertial(k,:)');
end
clear k;
errorDave = norm([errorMat_dave(1,:), errorMat_dave(2,:), errorMat_dave(3,:)])

%% Functions
function matrixTilda = skewSymmetric(vec)
% This function inputs a vector and returns a skew symmetrix matrix
matrixTilda = [0, -vec(3), vec(2);
               vec(3), 0, -vec(1);
               -vec(2), vec(1), 0;];
end