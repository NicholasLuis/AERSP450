clear; clc;

r1 = [1044.1; 7539.6; 1731.8];
r2 = [-229.6; 7389.3; 2491.0];
r3 = [-577.6; 7540.5; 1908.7]; 
r4 = [-3658.7; 6553.0; 2002.8];

r1hat = r1 ./ norm(r1);
r2hat = r2 ./ norm(r2);
r3hat = r3 ./ norm(r3);
r4hat = r4 ./ norm(r4);

% Not really h. Just the same direction as it
h = cross(r1hat, r3hat)

% Find which one does not equal zero
check1 = dot(h,r1hat)
check2 = dot(h,r2hat)
check3 = dot(h,r3hat)
check4 = dot(h,r4hat)