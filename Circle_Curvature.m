% ACCURACY ------------------------------------------------------------- %
% Determines amount of time steps performed and the amount of points used
% to create the shape(s).
A = 30; % Default 50
% ---------------------------------------------------------------------- %
% TIME PARAMETERS ------------------------------------------------------ %
T = 10;       % Length of time (default = 1)
S = 100;       % Amount of time steps
t = T/(S - 1); % Time interval
% ---------------------------------------------------------------------- %
% DISCRETE GRID -------------------------------------------------------- %
a = 0;             % Start of domain
b = 10;              % End of domain
n = A;               % Amount of points
x = linspace(a,b,n); % Set of equally spaced points
h = x(2) - x(1);     % Difference between points
[X,Y] = meshgrid(x); % Creates points to be used for phi function
% ---------------------------------------------------------------------- %
% PHI FUNCTION --------------------------------------------------------- %
%P = (0.5*(X-5).^2 + (4*(Y-5)).^2).^(0.5) - 3; % Circle
P = ((X-4).^2 + (Y-4).^2).^(0.5) - 2; % 2nd Circle
P = min(P,((X-6).^2 + (Y-6).^2).^(0.5) - 2); % 2nd Circle
% ---------------------------------------------------------------------- %
% INITIAL SPEED FUNCTION PARAMETER ------------------------------------- %
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
    shading interp                               % Interpolate for shade
    lighting gouraud                             % Gouraud light model
    light('Position',[0 0 100],'Style','local'); % Position local light
    pbaspect([1 1 1])
    grid on
    hold on
    patch([a b b a],[a a b b],[0 0 0 0],'k','FaceAlpha',0.5)
    hold off
    drawnow
  %  pause
    % UPDATE SPEED FUNCTION -------------------------------------------- %
    Pxx = diff(P([1 1:end end],:),2)/h^2;
    Pyy = diff(P(:,[1 1:end end])',2)'/h^2;
    Px = (P(3:end,:)-P(1:end-2,:))/(2*h); 
    Px = Px([1 1:end end],:);
    Py = (P(:,3:end)-P(:,1:end-2))/(2*h); 
    Py = Py(:,[1 1:end end]);
    Pxy = (Px(:,3:end)-Px(:,1:end-2))/(2*h); 
    Pxy = Pxy(:,[1 1:end end]);
    F = (Pxx.*Py.^2-2*Px.*Py.*Pxy+Pyy.*Px.^2)./(Px.^2+Py.^2).^1.5;
    F = min(max(F,-1/h),1/h);
    F = (1e-2)*F;
    F = -F;
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
    P = P - t*(max(F,0)*(gradPos) + min(F,0)*(gradNeg));
    % ------------------------------------------------------------------ %
end
% ---------------------------------------------------------------------- %