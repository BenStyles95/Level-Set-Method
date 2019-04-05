% TIME PARAMETERS ----------------------------------- %
T = 1;           % Length of time
S = 3;           % Amount of time steps
t = T / (S - 1); % Time interval
% --------------------------------------------------- %

% SHAPE ACCURACY ------------------------------------ %
n = 10000; % Amount of points to create shape 
% --------------------------------------------------- %

% SPEED FUNCTION ------------------------------------ %
% Code for functions of speed are placed here.
F = -1;
% --------------------------------------------------- %

% FUNCTION OF PHI ----------------------------------- %
% Code for functions of phi are placed here.
a = -50;
b = 50;
h = abs(b-a)/(n-1);
xValues = a:h:b;
phiValues = sin(xValues);
% --------------------------------------------------- %

% APPLY THE LEVEL SET METHOD'S UPWIND SCHEME -------- %

stepValues = zeros(n, S);    % Matrix to hold values
stepValues(:,1) = phiValues; % Store starting values

for i = 2:S
    diffPhi = diff(phiValues) / h;
    Dback = diffPhi([1 1:end]);
    Dforw = diffPhi([1:end end]);
    stepValues(:,i) = phiValues - t * F * (1+ max(Dforw,0).^2 + min(Dback,0).^2).^(1/2);
    phiValues = stepValues(:,i);
end

% PLOT STEPS ------------------------------------------------------------ %
% Plot the results on the same graph to track movement of the interface.
figure;
for i = 1:S
    plot(xValues, stepValues(:,i))
    hold on;
end
% ----------------------------------------------------------------------- %