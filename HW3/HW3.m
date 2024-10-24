% Made by Nicholas Luis (PSU ID 930841391)
% AERSP 450 HW 3
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
fprintf("Below are the values for perfect observations: \n")
v345 = getVelo(r3,r4,r5); norm(v345);
v246 = getVelo(r2,r4,r6); norm(v246);
v147 = getVelo(r1,r4,r7); norm(v147);
fprintf("Below is the δV (deviation from expected): \n")
dV345 = getDiff(v345 , V2true);
dV246 = getDiff(v246 , V2true);
dV147 = getDiff(v147 , V2true);
fprintf("Below is the %% Error (deviation from expected): \n")
pV345 = 100*dV345/norm(V2true);
pV246 = 100*dV246/norm(V2true);
pV147 = 100*dV147/norm(V2true);

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
fprintf("Below are the values for measured observations: \n")
v345 = getVelo(r3,r4,r5); norm(v345);
v246 = getVelo(r2,r4,r6); norm(v246);
v147 = getVelo(r1,r4,r7); norm(v147);
fprintf("Below is the δV (deviation from expected): \n")
dV345 = getDiff(v345 , V2true);
dV246 = getDiff(v246 , V2true);
dV147 = getDiff(v147 , V2true);
fprintf("Below is the %% Error (deviation from expected): \n")
pV345 = 100*dV345/norm(V2true);
pV246 = 100*dV246/norm(V2true);
pV147 = 100*dV147/norm(V2true);

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
fprintf("Below are the values for perfect observations: \n")
v345 = getVelo(r3,r4,r5); norm(v345);
v246 = getVelo(r2,r4,r6); norm(v246);
v147 = getVelo(r1,r4,r7); norm(v147);
fprintf("Below is the δV (deviation from expected): \n")
dV345 = getDiff(v345 , V2true);
dV246 = getDiff(v246 , V2true);
dV147 = getDiff(v147 , V2true);
fprintf("Below is the %% Error (deviation from expected): \n")
pV345 = 100*dV345/norm(V2true);
pV246 = 100*dV246/norm(V2true);
pV147 = 100*dV147/norm(V2true);

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
fprintf("Below are the values for measured observations: \n")
v345 = getVelo(r3,r4,r5); norm(v345);
v246 = getVelo(r2,r4,r6); norm(v246);
v147 = getVelo(r1,r4,r7); norm(v147);
fprintf("Below is the δV (deviation from expected): \n")
dV345 = getDiff(v345 , V2true);
dV246 = getDiff(v246 , V2true);
dV147 = getDiff(v147 , V2true);
fprintf("Below is the %% Error (deviation from expected): \n")
pV345 = 100*dV345/norm(V2true);
pV246 = 100*dV246/norm(V2true);
pV147 = 100*dV147/norm(V2true);

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
fprintf("Below are the values for perfect observations: \n")
v345 = getVelo(r3,r4,r5); norm(v345);
v246 = getVelo(r2,r4,r6); norm(v246);
v147 = getVelo(r1,r4,r7); norm(v147);
fprintf("Below is the δV (deviation from expected): \n")
dV345 = getDiff(v345 , V2true);
dV246 = getDiff(v246 , V2true);
dV147 = getDiff(v147 , V2true);
fprintf("Below is the %% Error (deviation from expected): \n")
pV345 = 100*dV345/norm(V2true);
pV246 = 100*dV246/norm(V2true);
pV147 = 100*dV147/norm(V2true);

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
fprintf("Below are the values for measured observations: \n")
v345 = getVelo(r3,r4,r5); norm(v345);
v246 = getVelo(r2,r4,r6); norm(v246);
v147 = getVelo(r1,r4,r7); norm(v147);
fprintf("Below is the δV (deviation from expected): \n")
dV345 = getDiff(v345 , V2true);
dV246 = getDiff(v246 , V2true);
dV147 = getDiff(v147 , V2true);
fprintf("Below is the %% Error (deviation from expected): \n")
pV345 = 100*dV345/norm(V2true);
pV246 = 100*dV246/norm(V2true);
pV147 = 100*dV147/norm(V2true);


%% Functions

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