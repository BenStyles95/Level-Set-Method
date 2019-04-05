% TIME PARAMETERS ------------------------------------------------------ %
T = 3.75;        % Length of time
S = 100;       % Amount of time steps
t = T/(S - 1); % Time interval
% ---------------------------------------------------------------------- %
% DOMAIN VALUES -------------------------------------------------------- %
a = 0;             % Start of domain 
b = 3;             % End of domain
n = 100;           % Number of points on domain
% BUILD DISCRETE GRID -------------------------------------------------- %
x = linspace(a,b,n); % Set of equally spaced points
h = x(2) - x(1);     % Difference between points
[X,Y] = meshgrid(x); % Creates points to be used for phi function    
% ---------------------------------------------------------------------- %
% LEVEL SET FUNCTION --------------------------------------------------- %
P = sin(pi*Y);
% ---------------------------------------------------------------------- %
% DYNAMIC INTERFACE TRACKING ------------------------------------------- %
figure('units','normalized','outerposition',[0 0 1 1]) % Full screen plot
for i = 1:S
    G = P;
    % GRAPH PLOTTING --------------------------------------------------- %
    set(gcf,'color','w')
    subplot(1,2,1)
    contourf(x,x,G,[0 0],'k-')
    set(gca,'fontsize',20)
    title(sprintf('Zero Level Set at t = %0.4f', (i-1)*t))
    axis([a b a b])
    pbaspect([1 1 1])
    subplot(1,2,2)
    surf(x,x,G,'EdgeAlpha',0.2)
    title('Level Set Function')
    set(gca,'fontsize',20)
    shading interp                         
    lighting gouraud                        
    light('Position',[0 0 max(G(:))+100],'Style','local');
    pbaspect([1 1 1])
    hold on
    patch([a b b a],[a a b b],[0 0 0 0],'k','FaceAlpha',0.5)
    hold off
    drawnow
    % ------------------------------------------------------------------ %
    % UPDATE SPEED FUNCTION -------------------------------------------- %
    F = -sin(pi*Y + i*t);
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