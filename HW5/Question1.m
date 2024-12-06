% AERSP 450 HW5 - Question 1 
% Made by Nicholas Luis (PSU ID: 930841391)

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

% Step 4
Star1 = [86.3143, 19.8000];
Star2 = [150.1887, 13.5660];
Star3 = [9.7843, 114.3137];
Star4 = [59.5000, 93.5000];
Star5 = [156.5000, 101.5000];