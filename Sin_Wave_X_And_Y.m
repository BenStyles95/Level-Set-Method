% ACCURACY ------------------------------------------------------------- %
A = 101; % Amount of points used to create front of interface
% ---------------------------------------------------------------------- %
% TIME PARAMETERS ------------------------------------------------------ %
T = 1;       % Length of time
S = 1000;    % Amount of time steps
t = T/(S-1); % Time interval
% ---------------------------------------------------------------------- %
% DISCRETE GRID -------------------------------------------------------- %
a = -2*pi;           % Start of domain
b = 2*pi;            % End of domain
n = A;               % Amount of points
x = linspace(a,b,n); % Set of equally spaced points
h = x(2)-x(1);       % Difference between points
[X,Y] = meshgrid(x); % Creates points to be used for phi function
% ---------------------------------------------------------------------- %
% PHI FUNCTION --------------------------------------------------------- %
P = sin(X);
% ---------------------------------------------------------------------- %
% DYNAMIC INTERFACE TRACKING ------------------------------------------- %
figure('units','normalized','outerposition',[0 0 1 1]) % Full screen plot
for i = 1:S
    time = i*t;
    % GRAPH PLOTTING --------------------------------------------------- %
    subplot(1,3,1)
    contourf(X,Y,P,[0 0],'k-')
    axis([a b a b])
    pbaspect([1 1 1])
    grid on
    subplot(1,3,2)
    contourf(X,Y,P)
    axis([a b a b])
    pbaspect([1 1 1])
    grid on
    subplot(1,3,3)
    surf(X,Y,P)
    % Lighting routine
    shading interp
    lighting gouraud                                       
    light('Position',[a a max(P(:))],'Style','local'); 
    pbaspect([1 1 1])
    grid on
    hold on
    patch([a b b a],[a a b b],[0 0 0 0],'k','FaceAlpha',0.5)
    hold off
    drawnow
    % UPDATE SPEED FUNCTION -------------------------------------------- %
    F = sin(Y+time);
    % LEVEL SET METHOD'S UPWIND SCHEME --------------------------------- %
    % Forward and backward differences between x points
    xDiff = diff(P)/h;
    xBackD = xDiff([1 1:end],:);
    xForwD = xDiff([1:end end],:);
    % Forward and backward differences between y points
    yDiff = diff(P')'/h;
    yBackD = yDiff(:,[1 1:end]);
    yForwD = yDiff(:,[1:end end]);
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