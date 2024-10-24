% Author: Melik Demirel
% mcd5703
% PSU ID: 952718091
% This code solves HW3 Problem 1
clc; close all; clear;

%% Givens

mu = 3.98600*10^14; % in SI [m]

% Case file names
caseFNS = {"Case1.mat", "Case2.mat", "Case3.mat"};
% Rtrue = True observations
% Rmeas = corrupted observations
% Tmeas = observation times  
% V2true = true velocity at the middle observation

%% Table generation

% Make a table and output it for each of the case files
for i = 1:length(caseFNS)
    generateCaseTable(mu, caseFNS{i}, i)
end

%% Functions

% Function to generate a table for each case file, given case i
function generateCaseTable(mu, file, i)
    fprintf('\n------------------\n')
    fprintf('----- CASE %i -----\n', i)
    fprintf('------------------\n')
    % Custom function for observation 3,4,5
    [trueV1, corrV1, delV1, pcErr1, dt1] = gibbsComp(mu, file, 3, 4, 5, i);
    % Custom function for observation 2,4,6
    [trueV2, corrV2, delV2, pcErr2, dt2] = gibbsComp(mu, file, 2, 4, 6, i);
    % Custom function for observation 1,4,7
    [trueV3, corrV3, delV3, pcErr3, dt3] = gibbsComp(mu, file, 1, 4, 7, i);

    % Align the data into a cell matrix
    data = {
    ' ',     'Est. Vel.','δv','% error',...
            'Est. Vel.','δv','% error','∆t1','∆t2';
    '3,4,5', v2s(trueV1), delV1(1), pcErr1(1), ...
                v2s(corrV1), delV1(2), pcErr1(2), dt1(1), dt1(2);
    '2,4,6', v2s(trueV2), delV2(1), pcErr2(1), ...
                v2s(corrV2), delV2(2), pcErr2(2), dt2(1), dt2(2);
    '1,4,7', v2s(trueV3), delV3(1), pcErr3(1), ...
                v2s(corrV3), delV3(2), pcErr3(2), dt3(1), dt3(2)
    };

    % Convert to table
    T = cell2table(data, 'VariableNames', ...
        {'Obs', 'Perfect Obs', ' ', '  ', 'Corrupted Obs', '   ', ...
                            '    ', '∆t1', '∆t2'});
    T
    filename = ['MelikGibbsTable.xlsx'];
    writetable(T, filename, 'Sheet', ['Case ', num2str(i)]);
end

% Velocity Vector to string
function out = v2s(vec)
    out = sprintf('[%.3f, %.3f, %.3f] km/s', ...
        vec(1)/1000, vec(2)/1000, vec(3)/1000);
end

% Function to extract the following data given a file
% names, mu, and o1, o2, and o3, and case i as simple integers 
% of the observation number:
% gt = gibbs true velocity (1x3 array)
% gc = gibbs corrupted velocity (1x3 array)
% err = error between the gibbs velocity and actual velocity
%       2 values for true and corrupted
% pc_error = percent error between the gibbs velocity and actual
%       2 values for true and corrupted
% dt = time difference between the observations
%       2 values for before and after time steo
function [gtv, gcv, err, pc_err, dt] = gibbsComp(mu, file, o1, o2, o3, i)
    % Load the data
    s = load(file); % in SI [m]
    % Use the custom gibbs function on the true data, 
    % then corrupted
    gtv = gibbs(mu,o1,o2,o3,s.Rtrue, ...
        [num2str(o1), ',', num2str(o2),',', num2str(o3), '-','True']);
    gcv = gibbs(mu,o1,o2,o3,s.RMeas, ...
        [num2str(o1), ',', num2str(o2),',', num2str(o3), '-','Corrupted']);
    % error calculation for true and corrupted
    err(1) = norm(gtv-s.V2true');
    err(2) = norm(gtv-s.V2true');
    pc_err(1) = norm(gtv-s.V2true') / norm(s.V2true) * 100;
    pc_err(2) = norm(gcv-s.V2true') / norm(s.V2true) * 100;
    % Delta t
    t = s.Tmeas;
    ti = t(o1);
    tm = t(o2);
    tf = t(o3);
    dti = tm - ti;
    dtf = tf - tm;
    dt = [dti, dtf];
end

% Given 3 observations, use the gibbs method on the data to find v2.
% mu is gravitation parameter in [m]
% o1, o2, and o3 are simple integers of the observation number
% in Data
% Data is a matrix of rows = number of observations and 3 columns
% for each vector direction
% CN is a string of the case name to be displayed
function v = gibbs(mu,o1,o2,o3,Data,CN)
    disp(CN);
    r1_ = Data(o1,:);
    r2_ = Data(o2,:);
    r3_ = Data(o3,:);

    % Gibbs Verification
    gibbsCheck = round(dot((r1_ ./ norm(r1_)), ...
        (cross(r2_,r3_)./norm(cross(r2_,r3_)))), 3);
    if gibbsCheck == 0
        disp('Gibbs Coplanar Verification Check Passed!')
    else
        error('Observation vectors are not coplanar!')
    end

    R1 = norm(r1_);
    R2 = norm(r2_);
    R3 = norm(r3_);

    % n, d, s calculation
    n = R1*cross(r2_,r3_) + R2*cross(r3_,r1_) + R3*cross(r1_,r2_);
    d = cross(r1_,r2_) + cross(r2_,r3_) + cross(r3_,r1_);
    s = r1_*(R2-R3) + r2_*(R3-R1) + r3_*(R1-R2);
    fprintf('\t n = [%.3e i, %.3e j, %.3e k] km.\n', ...
        n(1)/1000, n(2)/1000, n(3)/1000);
    fprintf('\t d = [%.3e i, %.3e j, %.3e k] km.\n', ...
        d(1)/1000, d(2)/1000, d(3)/1000);
    fprintf('\t s = [%.3e i, %.3e j, %.3e k] km.\n', ...
        s(1)/1000, s(2)/1000, s(3)/1000);

    % Gibbs calculation for v
    v = sqrt(mu / (norm(n)*norm(d)) ) * (cross(d,r2_)/R2 + s);
    fprintf('\t v = [%.3e i, %.3e j, %.3e k] km/s.\n', ...
        v(1)/1000, v(2)/1000, v(3)/1000);
end