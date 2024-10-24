% Made by Nicholas Luis (PSU ID 930841391)
% AERSP 450 HW 3 (Part B)

%% Case 1
% Getting the Orbital Elements (Uncomment to see the data)

[a, e, i, Omega, omega] = getOrbitalElements

% True Anomalies

%% Functions
function trueAnomaly = getF(E, R)
    % This function gets the true anomaly (in rad) given eccentricity and position vectors
    e = norm(E);
    r = norm(r);

    trueAnomaly = acos(dot(R,E) / (r*e)); 
end

function orbitalElements = getOrbitalElements(R,V)
    % This function gets the classical orbital elements given r and v vectors
    R = R / 1000; % Converting to km
    V = V / 1000; % Converting to km/s
    r = norm(R); 
    v = norm(V);
    MU = 3.986*(10^5); % km^3 / s^2

    % Getting semi major axis, a
    a = ((-2/MU)*(0.5*v^2 - MU/r))^(-1);
    
    % Intermediate Step: specific angular momentum
    h = cross(R,V);

    % Calculating eccentricity
    e = 

    orbitalElements = zeros(3,5); % 5 Columns, 3 Rows in each
    % Column 1 = a, Semi major axis (scalar stored in first element)
    % Column 2 = e, eccentricity (vector)
    % Column 3 = i, inclination (scalar stored in first element)
    % Column 4 = Ω, long. of the ascending node (scalar; first element)
    % Column 5 = ω, arg. of pariapsis (sclar; first element)
end