% Made by Nicholas Luis (PSU ID 930841391)
% AERSP 450 HW 3 (Part A)
clear; clc;

% Constants
MU = 3.986*(10^14); % m^3 / s^2

%% Case 1 - Perfect Observations
clear;
load('Case1.mat')

% Perfect Observations
r1 = Rtrue(1,:);
r2 = Rtrue(2,:);
r3 = Rtrue(3,:);
r4 = Rtrue(4,:);
r5 = Rtrue(5,:);
r6 = Rtrue(6,:);
r7 = Rtrue(7,:);
if coplanarCheck(r3,r4,r5)~=0, fprintf("r 3,4,5 are not coplanar!"), end
if coplanarCheck(r2,r4,r6)~=0, fprintf("r 2,4,6 are not coplanar!"), end
if coplanarCheck(r1,r4,r7)~=0, fprintf("r 1,4,7 are not coplanar!"), end

% Remove semicolons to print the data to copy + paste into excel
%fprintf("Below are the values for perfect observations: \n")
v345 = getVelo(r3,r4,r5); norm(v345);
v246 = getVelo(r2,r4,r6); norm(v246);
v147 = getVelo(r1,r4,r7); norm(v147);
%fprintf("Below is the δV (deviation from expected): \n")
dV345 = getDiff(v345 , V2true);
dV246 = getDiff(v246 , V2true);
dV147 = getDiff(v147 , V2true);
%fprintf("Below is the %% Error (deviation from expected): \n")
pV345 = 100*dV345/norm(V2true);
pV246 = 100*dV246/norm(V2true);
pV147 = 100*dV147/norm(V2true);

% PART B
% Getting Orbital Elements (Uncomment to output the data)

fprintf("Orbit elements for the 3-4-5 data set: \n")
orbitElements = getOrbitalElements(r4, v345);
a = orbitElements(1,1);
e = orbitElements(:,2);
e_mag = norm(e);
i = orbitElements(1, 3);
Omega = orbitElements(1, 4);
omega = orbitElements(1,5);

df = getF(e, r4) - getF(e, r3);
df = getF(e, r5) - getF(e, r4);

fprintf("Orbit elements for the 2-4-6 data set: \n")
orbitElements = getOrbitalElements(r4, v246);
a = orbitElements(1,1);
e = orbitElements(:,2);
e_mag = norm(e);
i = orbitElements(1, 3);
Omega = orbitElements(1, 4);
omega = orbitElements(1,5);

df = getF(e, r4) - getF(e, r2)
df = getF(e, r6) - getF(e, r4)

fprintf("Orbit elements for the 1-4-7 data set: \n")
orbitElements = getOrbitalElements(r4, v147);
a = orbitElements(1,1);
e = orbitElements(:,2);
e_mag = norm(e);
i = orbitElements(1, 3);
Omega = orbitElements(1, 4);
omega = orbitElements(1,5);

df = getF(e, r4) - getF(e, r1)
df = getF(e, r7) - getF(e, r4)

%% Case 1 - Corrupted Observations
clear;
load('Case1.mat')

% Corrupted Observations
r1 = RMeas(1,:);
r2 = RMeas(2,:);
r3 = RMeas(3,:);
r4 = RMeas(4,:);
r5 = RMeas(5,:);
r6 = RMeas(6,:);
r7 = RMeas(7,:);
if coplanarCheck(r3,r4,r5)~=0, fprintf("r 3,4,5 are not coplanar!"), end
if coplanarCheck(r2,r4,r6)~=0, fprintf("r 2,4,6 are not coplanar!"), end
if coplanarCheck(r1,r4,r7)~=0, fprintf("r 1,4,7 are not coplanar!"), end

% Remove semicolons to print the data to copy + paste into excel
%fprintf("Below are the values for measured observations: \n")
v345 = getVelo(r3,r4,r5); norm(v345);
v246 = getVelo(r2,r4,r6); norm(v246);
v147 = getVelo(r1,r4,r7); norm(v147);
%fprintf("Below is the δV (deviation from expected): \n")
dV345 = getDiff(v345 , V2true);
dV246 = getDiff(v246 , V2true);
dV147 = getDiff(v147 , V2true);
%fprintf("Below is the %% Error (deviation from expected): \n")
pV345 = 100*dV345/norm(V2true);
pV246 = 100*dV246/norm(V2true);
pV147 = 100*dV147/norm(V2true);

% PART B
% Getting Orbital Elements (Uncomment to output the data)

fprintf("Orbit elements for the 3-4-5 data set: \n")
orbitElements = getOrbitalElements(r4, v345);
a = orbitElements(1,1);
e = orbitElements(:,2);
e_mag = norm(e);
i = orbitElements(1, 3);
Omega = orbitElements(1, 4);
omega = orbitElements(1,5);

df = getF(e, r4) - getF(e, r3);
df = getF(e, r5) - getF(e, r4);

fprintf("Orbit elements for the 2-4-6 data set: \n")
orbitElements = getOrbitalElements(r4, v246);
a = orbitElements(1,1);
e = orbitElements(:,2);
e_mag = norm(e);
i = orbitElements(1, 3);
Omega = orbitElements(1, 4);
omega = orbitElements(1,5);

df = getF(e, r4) - getF(e, r2)
df = getF(e, r6) - getF(e, r4)

fprintf("Orbit elements for the 1-4-7 data set: \n")
orbitElements = getOrbitalElements(r4, v147);
a = orbitElements(1,1);
e = orbitElements(:,2);
e_mag = norm(e);
i = orbitElements(1, 3);
Omega = orbitElements(1, 4);
omega = orbitElements(1,5);

df = getF(e, r4) - getF(e, r1)
df = getF(e, r7) - getF(e, r4)

%% Case 2 - Perfect Observations
fprintf("\n\n\n-------BEGINNING OF CASE 2-------\n\n")
clear;
load('Case2.mat')

% Perfect Observations
r1 = Rtrue(1,:);
r2 = Rtrue(2,:);
r3 = Rtrue(3,:);
r4 = Rtrue(4,:);
r5 = Rtrue(5,:);
r6 = Rtrue(6,:);
r7 = Rtrue(7,:);
if coplanarCheck(r3,r4,r5)~=0, fprintf("r 3,4,5 are not coplanar!"), end
if coplanarCheck(r2,r4,r6)~=0, fprintf("r 2,4,6 are not coplanar!"), end
if coplanarCheck(r1,r4,r7)~=0, fprintf("r 1,4,7 are not coplanar!"), end

% Remove semicolons to print the data to copy + paste into excel
%fprintf("Below are the values for perfect observations: \n")
v345 = getVelo(r3,r4,r5); norm(v345);
v246 = getVelo(r2,r4,r6); norm(v246);
v147 = getVelo(r1,r4,r7); norm(v147);
%fprintf("Below is the δV (deviation from expected): \n")
dV345 = getDiff(v345 , V2true);
dV246 = getDiff(v246 , V2true);
dV147 = getDiff(v147 , V2true);
%fprintf("Below is the %% Error (deviation from expected): \n")
pV345 = 100*dV345/norm(V2true);
pV246 = 100*dV246/norm(V2true);
pV147 = 100*dV147/norm(V2true);

% PART B
% Getting Orbital Elements (Uncomment to output the data)

fprintf("Orbit elements for the 3-4-5 data set: \n")
orbitElements = getOrbitalElements(r4, v345);
a = orbitElements(1,1);
e = orbitElements(:,2);
e_mag = norm(e);
i = orbitElements(1, 3);
Omega = orbitElements(1, 4);
omega = orbitElements(1,5);

df = getF(e, r4) - getF(e, r3);
df = getF(e, r5) - getF(e, r4);

fprintf("Orbit elements for the 2-4-6 data set: \n")
orbitElements = getOrbitalElements(r4, v246);
a = orbitElements(1,1);
e = orbitElements(:,2);
e_mag = norm(e);
i = orbitElements(1, 3);
Omega = orbitElements(1, 4);
omega = orbitElements(1,5);

df = getF(e, r4) - getF(e, r2)
df = getF(e, r6) - getF(e, r4)

fprintf("Orbit elements for the 1-4-7 data set: \n")
orbitElements = getOrbitalElements(r4, v147);
a = orbitElements(1,1);
e = orbitElements(:,2);
e_mag = norm(e);
i = orbitElements(1, 3);
Omega = orbitElements(1, 4);
omega = orbitElements(1,5);

df = getF(e, r4) - getF(e, r1)
df = getF(e, r7) - getF(e, r4)

%% Case 2 - Corrupted Observations
clear;
load('Case2.mat')

% Corrupted Observations
r1 = RMeas(1,:);
r2 = RMeas(2,:);
r3 = RMeas(3,:);
r4 = RMeas(4,:);
r5 = RMeas(5,:);
r6 = RMeas(6,:);
r7 = RMeas(7,:);
if coplanarCheck(r3,r4,r5)~=0, fprintf("r 3,4,5 are not coplanar!"), end
if coplanarCheck(r2,r4,r6)~=0, fprintf("r 2,4,6 are not coplanar!"), end
if coplanarCheck(r1,r4,r7)~=0, fprintf("r 1,4,7 are not coplanar!"), end

% Remove semicolons to print the data to copy + paste into excel
%fprintf("Below are the values for measured observations: \n")
v345 = getVelo(r3,r4,r5); norm(v345);
v246 = getVelo(r2,r4,r6); norm(v246);
v147 = getVelo(r1,r4,r7); norm(v147);
%fprintf("Below is the δV (deviation from expected): \n")
dV345 = getDiff(v345 , V2true);
dV246 = getDiff(v246 , V2true);
dV147 = getDiff(v147 , V2true);
%fprintf("Below is the %% Error (deviation from expected): \n")
pV345 = 100*dV345/norm(V2true);
pV246 = 100*dV246/norm(V2true);
pV147 = 100*dV147/norm(V2true);

% PART B
% Getting Orbital Elements (Uncomment to output the data)

fprintf("Orbit elements for the 3-4-5 data set: \n")
orbitElements = getOrbitalElements(r4, v345);
a = orbitElements(1,1);
e = orbitElements(:,2);
e_mag = norm(e);
i = orbitElements(1, 3);
Omega = orbitElements(1, 4);
omega = orbitElements(1,5);

df = getF(e, r4) - getF(e, r3);
df = getF(e, r5) - getF(e, r4);

fprintf("Orbit elements for the 2-4-6 data set: \n")
orbitElements = getOrbitalElements(r4, v246);
a = orbitElements(1,1);
e = orbitElements(:,2);
e_mag = norm(e);
i = orbitElements(1, 3);
Omega = orbitElements(1, 4);
omega = orbitElements(1,5);

df = getF(e, r4) - getF(e, r2)
df = getF(e, r6) - getF(e, r4)

fprintf("Orbit elements for the 1-4-7 data set: \n")
orbitElements = getOrbitalElements(r4, v147);
a = orbitElements(1,1);
e = orbitElements(:,2);
e_mag = norm(e);
i = orbitElements(1, 3);
Omega = orbitElements(1, 4);
omega = orbitElements(1,5);

df = getF(e, r4) - getF(e, r1)
df = getF(e, r7) - getF(e, r4)

%% Case 3 - Perfect Observations
fprintf("\n\n\n-------BEGINNING OF CASE 3-------\n\n")
clear;
load('Case3.mat')

% Perfect Observations
r1 = Rtrue(1,:);
r2 = Rtrue(2,:);
r3 = Rtrue(3,:);
r4 = Rtrue(4,:);
r5 = Rtrue(5,:);
r6 = Rtrue(6,:);
r7 = Rtrue(7,:);
if coplanarCheck(r3,r4,r5)~=0, fprintf("r 3,4,5 are not coplanar!"), end
if coplanarCheck(r2,r4,r6)~=0, fprintf("r 2,4,6 are not coplanar!"), end
if coplanarCheck(r1,r4,r7)~=0, fprintf("r 1,4,7 are not coplanar!"), end

% Remove semicolons to print the data to copy + paste into excel
%fprintf("Below are the values for perfect observations: \n")
v345 = getVelo(r3,r4,r5); norm(v345);
v246 = getVelo(r2,r4,r6); norm(v246);
v147 = getVelo(r1,r4,r7); norm(v147);
%fprintf("Below is the δV (deviation from expected): \n")
dV345 = getDiff(v345 , V2true);
dV246 = getDiff(v246 , V2true);
dV147 = getDiff(v147 , V2true);
%fprintf("Below is the %% Error (deviation from expected): \n")
pV345 = 100*dV345/norm(V2true);
pV246 = 100*dV246/norm(V2true);
pV147 = 100*dV147/norm(V2true);

% PART B
% Getting Orbital Elements (Uncomment to output the data)

fprintf("Orbit elements for the 3-4-5 data set: \n")
orbitElements = getOrbitalElements(r4, v345);
a = orbitElements(1,1);
e = orbitElements(:,2);
e_mag = norm(e);
i = orbitElements(1, 3);
Omega = orbitElements(1, 4);
omega = orbitElements(1,5);

df = getF(e, r4) - getF(e, r3);
df = getF(e, r5) - getF(e, r4);

fprintf("Orbit elements for the 2-4-6 data set: \n")
orbitElements = getOrbitalElements(r4, v246);
a = orbitElements(1,1);
e = orbitElements(:,2);
e_mag = norm(e);
i = orbitElements(1, 3);
Omega = orbitElements(1, 4);
omega = orbitElements(1,5);

df = getF(e, r4) - getF(e, r2)
df = getF(e, r6) - getF(e, r4)

fprintf("Orbit elements for the 1-4-7 data set: \n")
orbitElements = getOrbitalElements(r4, v147);
a = orbitElements(1,1);
e = orbitElements(:,2);
e_mag = norm(e);
i = orbitElements(1, 3);
Omega = orbitElements(1, 4);
omega = orbitElements(1,5);

df = getF(e, r4) - getF(e, r1)
df = getF(e, r7) - getF(e, r4)

%% Case 3 - Corrupted Observations
load('Case3.mat')

% Corrupted Observations
r1 = RMeas(1,:);
r2 = RMeas(2,:);
r3 = RMeas(3,:);
r4 = RMeas(4,:);
r5 = RMeas(5,:);
r6 = RMeas(6,:);
r7 = RMeas(7,:);
if coplanarCheck(r3,r4,r5)~=0, fprintf("r 3,4,5 are not coplanar!"), end
if coplanarCheck(r2,r4,r6)~=0, fprintf("r 2,4,6 are not coplanar!"), end
if coplanarCheck(r1,r4,r7)~=0, fprintf("r 1,4,7 are not coplanar!"), end

% Remove semicolons to print the data to copy + paste into excel
%fprintf("Below are the values for measured observations: \n")
v345 = getVelo(r3,r4,r5); norm(v345);
v246 = getVelo(r2,r4,r6); norm(v246);
v147 = getVelo(r1,r4,r7); norm(v147);
%fprintf("Below is the δV (deviation from expected): \n")
dV345 = getDiff(v345 , V2true);
dV246 = getDiff(v246 , V2true);
dV147 = getDiff(v147 , V2true);
%fprintf("Below is the %% Error (deviation from expected): \n")
pV345 = 100*dV345/norm(V2true);
pV246 = 100*dV246/norm(V2true);
pV147 = 100*dV147/norm(V2true);

% PART B
% Getting Orbital Elements (Uncomment to output the data)

fprintf("Orbit elements for the 3-4-5 data set: \n")
orbitElements = getOrbitalElements(r4, v345);
a = orbitElements(1,1);
e = orbitElements(:,2);
e_mag = norm(e);
i = orbitElements(1, 3);
Omega = orbitElements(1, 4);
omega = orbitElements(1,5);

df = getF(e, r4) - getF(e, r3);
df = getF(e, r5) - getF(e, r4);

fprintf("Orbit elements for the 2-4-6 data set: \n")
orbitElements = getOrbitalElements(r4, v246);
a = orbitElements(1,1);
e = orbitElements(:,2);
e_mag = norm(e);
i = orbitElements(1, 3);
Omega = orbitElements(1, 4);
omega = orbitElements(1,5);

df = getF(e, r4) - getF(e, r2)
df = getF(e, r6) - getF(e, r4)

fprintf("Orbit elements for the 1-4-7 data set: \n")
orbitElements = getOrbitalElements(r4, v147);
a = orbitElements(1,1);
e = orbitElements(:,2);
e_mag = norm(e);
i = orbitElements(1, 3);
Omega = orbitElements(1, 4);
omega = orbitElements(1,5);

df = getF(e, r4) - getF(e, r1)
df = getF(e, r7) - getF(e, r4)

%% Functions
function trueAnomaly = getF(E, R)
    % This function gets the true anomaly (in deg) given eccentricity and position vectors
    e = norm(E);
    r = norm(R);

    trueAnomaly = acos(dot(R,E) / (r*e)) * (180/pi); 
end

function orbitalElements = getOrbitalElements(R,V)
    % This function gets the classical orbital elements given r and v vectors
    R = R / 1000; % Converting to km
    V = V / 1000; % Converting to km/s
    r = norm(R); 
    v = norm(V);
    MU = 3.986*(10^5); % km^3 / s^2
    i = [1 0 0]; j = [0 1 0]; k = [0 0 1];

    % Getting semi major axis, a
    a = ((-2/MU)*(0.5*v^2 - MU/r))^(-1);
    
    % Intermediate Steps
    h = cross(R,V);
    n = cross(k,h);

    % Calculating eccentricity
    e = (cross(V,h)/MU) - (R/r);

    % Calculating inclination
    I = acos(dot(k,h/norm(h))) * (180/pi);

    % Calculating longitude of ascending node
    Omega = atan2(dot(j,n/norm(n)), dot(j,n/norm(n))) * (180/pi);

    % Calculating argument of periapsis
    omega = acos(dot(n/norm(n),e/norm(e))) * (180/pi);

    orbitalElements = NaN(3,5); % 5 Columns, 3 Rows in each
    % Column 1 = a, Semi major axis (scalar stored in first element)
    orbitalElements(1,1) = a;
    % Column 2 = e, eccentricity (vector)
    orbitalElements(:,2) = e';
    % Column 3 = i, inclination (scalar stored in first element)
    orbitalElements(1,3) = I;
    % Column 4 = Ω, long. of the ascending node (scalar; first element)
    orbitalElements(1,4) = Omega;
    % Column 5 = ω, arg. of pariapsis (sclar; first element)
    orbitalElements(1,5) = omega;
end

function delta = getDiff(va, vb)
    % This function gets magnitude difference between two vectors of length 3
    delta = sqrt( (va(1)-vb(1))^2 + (va(2)-vb(2))^2 + (va(3)-vb(3))^2 );
end

function vOut = getVelo(ra, rb, rc)
    MU = 3.986*(10^14); % m^3 / s^2 
    % This function gets velocity at the middle position
    n = (norm(ra).*cross(rb,rc)) + (norm(rb).*cross(rc,ra)) + (norm(rc).*cross(ra,rb));
    d = cross(ra,rb) + cross(rb,rc) + cross(rc,ra);
    s = (ra.*(norm(rb)-norm(rc))) + (rb.*(norm(rc)-norm(ra))) + (rc.*(norm(ra)-norm(rb)));

    vOut = ((cross(d,rb)/norm(rb))+s).*sqrt(MU/(norm(n)*norm(d)));
end

function boolean = coplanarCheck(ra, rb, rc)
    % This function checks if all the vectors are coplanar (within reason)

    % Converting to unit vectors
    ra = ra ./ norm(ra);
    rb = rb ./ norm(rb);
    rc = rc ./ norm(rc);

    boolean =  round( dot(ra, cross(rb,rc)), 10 ); % Rounds if close enough to 0
end