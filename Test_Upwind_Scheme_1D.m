% TIME PARAMETERS ----------------------------------- %
T = 1;           % Length of time
S = 100;         % Amount of time steps
t = T / (S - 1); % Time interval
% --------------------------------------------------- %
% SHAPE ACCURACY ------------------------------------ %
n = 100; % Amount of points to create shape 
% --------------------------------------------------- %
% SPEED FUNCTION ------------------------------------ %
F = -1;
% --------------------------------------------------- %
% FUNCTION OF PHI ----------------------------------- %
% Code for functions of phi are placed here.
a = 0;
b = 1;
h = abs(b - a) / (n - 1);
xValues = a:h:b;
yValues = zeros;
for i = 1:n
    if xValues(i) <= 1/2
        yValues(i) = 1/2 - xValues(i);
    elseif xValues(i) > 1/2
        yValues(i) = xValues(i) - 1/2;
    end
end
% --------------------------------------------------- %
% APPLY THE LEVEL SET METHOD'S UPWIND SCHEME -------- %
stepMatrix = zeros(n, S);    % Matrix to hold values
stepMatrix(:,1) = yValues;   % Store starting values
for i = 2:S
    diffY = diff(yValues) / h;
    dXBack = diffY([1 1:end]);
    dXForw = diffY([1:end end]);
    stepMatrix(:,i) = yValues - t * F * (min(dXBack,0).^2 + max(dXForw,0).^2).^(1/2);
    yValues = stepMatrix(:,i);
end
% PLOT STEPS ------------------------------------------------------------ %
% Plot the results on the same graph to track movement of the interface.
figure;
for i = 1:S
    plot(xValues, stepMatrix(:,i))
    axis([0 1 0 1.5]) 
    drawnow
end
% ----------------------------------------------------------------------- %