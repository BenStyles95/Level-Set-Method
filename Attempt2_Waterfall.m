% ACCURACY ------------------------------------------------------------- %
% Determines amount of time steps performed and the amount of points used
% to create the shape(s). Lower values have better performance, higher
% values have better accuracy.
A = 50; % Default 50
% ---------------------------------------------------------------------- %
% TIME PARAMETERS ------------------------------------------------------ %
T = 5;         % Length of time (default = 1)
S = 100*T;     % Amount of time steps
t = T/(S - 1); % Time interval
% ---------------------------------------------------------------------- %
% DISCRETE GRID -------------------------------------------------------- %
a = -10;             % Start of domain
b = 10;              % End of domain
n = A;               % Amount of points
x = linspace(a,b,n); % Set of equally spaced points
y = x;
h = x(2) - x(1);     % Difference between points
[X,Y] = meshgrid(x,y); % Creates points to be used for phi function
% ---------------------------------------------------------------------- %
% PHI FUNCTION --------------------------------------------------------- %
P = zeros(n,n); 
arrayPi = linspace(0,2*pi,n);
for i = 1:n
    P(:,i) = sin(arrayPi);
end
% ---------------------------------------------------------------------- %
% DYNAMIC INTERFACE TRACKING ------------------------------------------- %
figure('units','normalized','outerposition',[0 0 1 1]) % Full screen plot
for i = 1:S
    % GRAPH PLOTTING --------------------------------------------------- %
    subplot(1,3,1)
    contourf(x,y,P,[0 0],'k-')
    axis([a b a b])
    pbaspect([1 1 1])
    grid on
    subplot(1,3,2)
    contourf(x,y,P)
    axis([a b a b])
    pbaspect([1 1 1])
    grid on
    subplot(1,3,3)
    surf(x,y,P,'EdgeAlpha',0.2)
    shading interp                               % Interpolate for shade
    lighting gouraud                             % Gouraud light model
    light('Position',[0 0 200],'Style','local'); % Position local light
    pbaspect([1 1 1])
    grid on  
    hold on
    patch([a b b a],[a a b b],[0 0 0 0],'k','FaceAlpha',0.5)
    hold off
    drawnow
    % UPDATE SPEED FUNCTION -------------------------------------------- %
    for j = 1:n
        F(:,j) = sin(arrayPi+i*t*10);
    end
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