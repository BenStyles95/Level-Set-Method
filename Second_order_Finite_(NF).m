% ACCURACY ------------------------------------------------------------- %
% Determines amount of time steps performed and the amount of points used
% to create the shape(s).
A = 50; % Default 50
% ---------------------------------------------------------------------- %
% TIME PARAMETERS ------------------------------------------------------ %
T = 1;         % Length of time (default = 1)
S = T*A;       % Amount of time steps
t = T/(S - 1); % Time interval
% ---------------------------------------------------------------------- %
% DISCRETE GRID -------------------------------------------------------- %
a = -10;             % Start of domain
b = 10;              % End of domain
n = A;               % Amount of points
x = linspace(a,b,n); % Set of equally spaced points
h = x(2) - x(1);     % Difference between points
[X,Y] = meshgrid(x); % Creates points to be used for phi function
% ---------------------------------------------------------------------- %
% PHI FUNCTION --------------------------------------------------------- %
P = ((X).^2 + (Y).^2).^(0.5) - 1; % Circle
% ---------------------------------------------------------------------- %
% INITIAL SPEED FUNCTION PARAMETER ------------------------------------- %
F = 1;
% ---------------------------------------------------------------------- %
% DYNAMIC INTERFACE TRACKING ------------------------------------------- %
figure('units','normalized','outerposition',[0 0 1 1]) % Full screen plot
for i = 1:S
    % GRAPH PLOTTING --------------------------------------------------- %
    subplot(1,3,1)
    contourf(x,x,-P,[0 0],'k-')
    axis([a b a b])
    pbaspect([1 1 1])
    grid on
    subplot(1,3,2)
    contourf(x,x,-P)
    axis([a b a b])
    pbaspect([1 1 1])
    grid on
    subplot(1,3,3)
    surf(x,x,-P,'EdgeAlpha',0.2)
    pbaspect([1 1 1])
    grid on
    hold on
    patch([a b b a],[a a b b],[0 0 0 0],'k','FaceAlpha',0.5)
    hold off
    drawnow
    % UPDATE SPEED FUNCTION -------------------------------------------- %
    F = F;
    % LEVEL SET METHOD'S UPWIND SCHEME --------------------------------- %
    % Forward and backward differences between x points
    xDiff = diff(P)/h;
    xBackD = xDiff([1 1:end],:);
    xForwD = xDiff([1:end end],:);
    
    xDiff2nd = diff(P,2);                       % 2nd order finite difference
    xBackD2nd = xDiff2nd([1 1 1:end],:);        % Backward 2nd order finite difference
    xForwD2nd = xDiff2nd([1:end end end],:);    % Forward 2nd order finite difference
    
    xBackD = xBackD + h*xBackD2nd;
    xForwD = xForwD - h*xForwD2nd;
    
    xDiff1stPoint = (P(1,:) - P(end,:))/h;
    
    xBackD(1,:) = xDiff1stPoint;
    xForwD(A,:) = xDiff1stPoint;
    % Forward and backward differences between y points
    yDiff = diff(P')'/h;
    yBackD = yDiff(:,[1 1:end]);
    yForwD = yDiff(:,[1:end end]);
    
    yDiff2nd = diff(P',2)'/h;                   %#ok<*UDIM> % 2nd order finite difference
    yBackD2nd = yDiff2nd(:,[1 1 1:end]);        % Backward 2nd order finite difference
    yForwD2nd = yDiff2nd(:,[1:end end end]);    % Forward 2nd order finite difference
    
    yBackD = yBackD + h*yBackD2nd;
    yForwD = yForwD - h*yForwD2nd;
    
    yDiff1stPoint = (P(:,1)' - P(:,end)')'/h;
    
    xBackD(:,1) = yDiff1stPoint;
    xForwD(:,A) = yDiff1stPoint;
    % Select the appropriate finite difference method (upwind scheme)
    gradPos = (max(xBackD,0).^2 + min(xForwD,0).^2 + ...
               max(yBackD,0).^2 + min(yForwD,0).^2).^(1/2);
    gradNeg = (min(xBackD,0).^2 + max(xForwD,0).^2 + ...
               min(yBackD,0).^2 + max(yForwD,0).^2).^(1/2);
    % Update points using the level set equation       
    P = P - t*(max(F,0)*gradPos + min(F,0)*gradNeg);
    % ------------------------------------------------------------------ %
end
% ---------------------------------------------------------------------- %