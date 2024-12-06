% AERSP 450 HW5 - Question 1 
% Made by Nicholas Luis (PSU ID: 930841391)

clc; clear; close all;

%% Part A

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

%% Part B
StarI1 = [-0.921069884293268 -0.342599924017704 0.185082577005643];
StarI2= [-0.476980639452282 -0.711962047538100 0.515363476056510];
StarI3= [-0.496592767571065 0.816782636490539 -0.293703503424244];
StarI4= [-0.736236774114155 0.676644212846351 0.010393347079969];
StarI5 = [ 0.304668730748568 -0.299446239108458 0.904161995655567];
starCoordsInertial = [StarI1; StarI2; StarI3; StarI4; StarI5]