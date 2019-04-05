% TIME PARAMETERS ------------------------------------------------------ %
T = 1;         % Length of time
S = 25*T;      % Amount of time steps
t = T/(S - 1); % Time interval
% ---------------------------------------------------------------------- %
% DISCRETE GRID -------------------------------------------------------- %
a = -10;               % Start of domain 
b = 10;                % End of domain
n = S / T;               % Amount of points
x = linspace(a,b,n);   % Set of equally spaced points
h = x(2) - x(1);       % Difference between points
[X,Y] = meshgrid(x);   % Creates points to be used for    
% ---------------------------------------------------------------------- %
% CO-ORDINATES OF SHAPE ------------------------------------------------ %
P = zeros(n,n,S);
P(:,:,1) = ((X).^2 + (Y).^2).^(1/2) - 1; % Circle at (0,0) with r = 1
% ---------------------------------------------------------------------- %
% SPEED FUNCTION ------------------------------------------------------- %
F = 1;
% ---------------------------------------------------------------------- %
% LEVEL SET METHOD'S UPWIND SCHEME ------------------------------------- %
for i = 1:(S - 1)    
    % Forward and backward differences between x points
    xDiff = diff(P(:,:,i))/h;
    xBackD = xDiff([1 1:end], :);
    xForwD = xDiff([1:end end], :);
    % Forward and backward differences between y points
    yDiff = diff(P(:,:,i)')'/h;
    yBackD = yDiff(:, [1 1:end]);
    yForwD = yDiff(:, [1:end end]); 
    % Select the appropriate finite difference method (upwind scheme)
    gradPos = (max(xBackD,0).^2 + min(xForwD,0).^2 + ...
               max(yBackD,0).^2 + min(yForwD,0).^2).^(1/2);
    gradNeg = (min(xBackD,0).^2 + max(xForwD,0).^2 + ...
               min(yBackD,0).^2 + max(yForwD,0).^2).^(1/2);           
    % Update points using the level set equation
    P(:,:,i+1) = P(:,:,i) - t*(max(F,0)*gradPos + min(F,0)*gradNeg);
end
% ---------------------------------------------------------------------- %
% GRAPH PLOTTING ------------------------------------------------------- %
figure('units','normalized','outerposition',[0 0 1 1]) % Full screen plot
for i = 1:S
    subplot(1,3,1) % 1st plot tracks shape's topology
    contourf(x,x,-(P(:,:,i)),[0 0],'k-')
    axis([a b a b])
    pbaspect([1 1 1])
    grid on
    subplot(1,3,2) % 2nd plot tracks contours of the shape
    contourf(x,x,-(P(:,:,i)))
    axis([a b a b])
    pbaspect([1 1 1])
    grid on
    subplot(1,3,3) % 3rd plot tracks shape's growth against time
    surf(x, x, -(P(:,:,i)), 'EdgeAlpha', 0.2)
    pbaspect([1 1 1])
    grid on
    hold on
    patch([a b b a],[a a b b],[0 0 0 0],'k','FaceAlpha',0.5)
    hold off
    drawnow
end
% ---------------------------------------------------------------------- %