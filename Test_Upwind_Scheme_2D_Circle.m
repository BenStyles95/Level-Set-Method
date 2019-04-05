% TIME PARAMETERS ------------------------------------------------------ %
T = 1;           % Length of time
S = 100;         % Amount of time steps
t = T / (S - 1); % Time interval
% ---------------------------------------------------------------------- %
% SHAPE ACCURACY ------------------------------------------------------- %
n = 100; % Amount of points to create shape 
% ---------------------------------------------------------------------- %
% SPEED FUNCTION ------------------------------------------------------- %
F = -1;
% ---------------------------------------------------------------------- %
% CO-ORDINATES OF SHAPE ------------------------------------------------ %
a = 0; b = 2*pi; h = abs(b-a)/(n-1);
ang = a:h:b;
xC = 3; yC = 3; r = 3;
xValues = r*cos(ang) + xC;
yValues = r*sin(ang) + yC;
% ---------------------------------------------------------------------- %
% APPLY THE LEVEL SET METHOD'S UPWIND SCHEME --------------------------- %
stepMatrixX = zeros(n, S+1);    % Matrix to hold values
stepMatrixX(:,1) = xValues;    % Store starting values
stepMatrixY = zeros(n, S+1);
stepMatrixY(:,1) = yValues;
for i = 1:S
    diffX = diff(stepMatrixX(:,i)) / h;
    dXBack = diffX([1 1:end]);
    dXForw = diffX([1:end end]);
    diffY = diff(stepMatrixY(:,i)) / h;
    dYBack = diffY([1 1:end]);
    dYForw = diffY([1:end end]);
    
    GradPlus=(max(dXBack,0).^2 + min(dXForw,0).^2 + max(dYBack,0).^2 + min(dYForw,0).^2).^(1/2);
    GradMinus=(min(dXBack,0).^2 + max(dXForw,0).^2 + min(dYBack,0).^2 + max(dYForw,0).^2).^(1/2);
    
    stepMatrixX(:,i+1) = stepMatrixX(:,i) - t * (max(F,0)*GradPlus+min(F,0)*GradMinus);
    stepMatrixY(:,i+1) = stepMatrixY(:,i)- t * (max(F,0)*GradPlus+min(F,0)*GradMinus); 
end
% ---------------------------------------------------------------------- %
% PLOT STEPS ----------------------------------------------------------- %
% Plot the results on the same graph to track movement of the interface.
figure;
for i = 1:S
    plot(stepMatrixX(:,i), stepMatrixY(:,i))
    axis([a-5 b+5 a-5 b+5])
    drawnow
end
% ---------------------------------------------------------------------- %