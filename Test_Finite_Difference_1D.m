% USER INPUT
L = 26;     % Amount of levels
n = 101;    % Amount of points      

% SPEED FUNCTION
% This is where a given speed function would be coded into the script. For 
% this case, I have made speed function equal to 1.
F = -1;

% FUNCTION FOR PHI
% For this case, phi(x) = {1/2 - x for x <= 1/2, x - 1/2 for x > 1/2} over
% the domain 0 < x < 1.
a = 0;
b = 1;
h = abs(b-a)/(n-1);
T = 0.25;
t = T/(L - 1);  % Time interval

xNums = a:h:b;
phiValues = zeros;
for i = 1:n
    if xNums(i) <= 1/2
        phiValues(i) = 1/2 - xNums(i);
    elseif xNums(i) > 1/2
        phiValues(i) = xNums(i) - 1/2;
    end
end

% APPLY LEVEL SET FORMULA
% For this case, the use of finite difference approximations were used.
phiNewValues = zeros(n, L);     % Matrix to hold phi values at each level
phiNewValues(:,1) = phiValues;  % Store starting values of phi (0 level)
for j = 2:L
    % Finite difference approximations to find n+1 phi values 
    phiNewValues(1, j) = phiValues(1) - t * F * (1 + ((phiValues(2) - phiValues(1)) / (h) ).^2 ).^(1/2);           % Forward difference
    for i = 2:n-1
        phiNewValues(i, j) = phiValues(i) - t * F * (1 + ((phiValues(i+1) - phiValues(i-1)) / (2*h) ).^2 ).^(1/2); % Central difference
    end
    phiNewValues(n, j) = phiValues(n) - t * F * (1 + ((phiValues(n) - phiValues(n-1)) / (h)).^2).^(1/2);           % Backward difference
    % Stores new phi values as starting values for next iteration
    phiValues = phiNewValues(:,j);
end

% CREATE ARRAY FOR ALL TIME STEPS
% Finding t(n) for all n.
tValues = zeros;    % Array to hold t(n) values
for i = 1:n
    tValues(i) = (i-1)*t;
end

% PLOT LEVELS
% Plot the results on the same graph to track movement of the interface.
figure;
for i = 1:L
  plot(xNums, phiNewValues(:,i))
  axis([0 1 0 1.25])  
  drawnow
end