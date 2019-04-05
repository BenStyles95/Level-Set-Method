% ACCURACY ------------------------------------------------------------- %
% Determines amount of time steps performed and the amount of points used 
% to create the shape(s).
A = 100;
% ---------------------------------------------------------------------- %
% TIME PARAMETERS ------------------------------------------------------ %
T = 1;         % Length of time
S = 100;       % Amount of time steps
t = T/(S - 1); % Time interval
% ---------------------------------------------------------------------- %
% DISCRETE GRID -------------------------------------------------------- %
a = -5;             % Start of domain 
b = 5;              % End of domain
n = A;               % Amount of points
x = linspace(a,b,n); % Set of equally spaced points
h = x(2) - x(1);     % Difference between points
[X,Y] = meshgrid(x); % Creates points to be used for phi function    
% ---------------------------------------------------------------------- %
% PHI FUNCTION --------------------------------------------------------- %%
c = 4/3; % Centre
r = 1; % Radius
P = ((X + c).^2 + (Y + c).^2).^(1/2) - r;       % 1st circle
P = min(P, ((X - c).^2 + (Y - c).^2).^(1/2) - r); % 2nd circle
P = min(P, ((X + c).^2 + (Y - c).^2).^(1/2) - r); % 3rd circle
P = min(P, ((X - c).^2 + (Y + c).^2).^(1/2) - r); % 4th circle
% ---------------------------------------------------------------------- %
% SPEED FUNCTION ------------------------------------------------------- %
F = 1;
% ---------------------------------------------------------------------- %
% DYNAMIC INTERFACE TRACKING ------------------------------------------- %
figure('units','normalized','outerposition',[0 0 1 1]) % Full screen plot
for i = 1:S
    % GRAPH PLOTTING --------------------------------------------------- %
    set(gcf,'color','w')
    subplot(1,2,1)
    contourf(x,x,-P,[0 0],'k-')
    set(gca,'fontsize',16)
    title(sprintf('Front at t = %0.4f', (i-1)*t))
    axis([a b a b])
    pbaspect([1 1 1])
    subplot(1,2,2)
    surf(x,x,-P,'EdgeAlpha',0.2)
    title('Surface Plot of Level Set Function')
    set(gca,'fontsize',16)
    shading interp                               % Interpolate for shade
    lighting gouraud                             % Gouraud light model
    light('Position',[0 0 max(P(:))+100],'Style','local'); % Position local light
    pbaspect([1 1 1])
    hold on
    patch([a b b a],[a a b b],[0 0 0 0],'k','FaceAlpha',0.5)
    hold off
    drawnow
    %pause
    % ------------------------------------------------------------------ %
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