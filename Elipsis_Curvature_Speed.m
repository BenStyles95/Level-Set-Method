% TIME PARAMETERS ------------------------------------------------------ %
T = 3;         % Length of time
S = 300;      % Amount of time steps
t = T/(S - 1); % Time interval
% ---------------------------------------------------------------------- %
% DOMAIN VALUES -------------------------------------------------------- %
a = -5;            % Start of domain 
b = 5;             % End of domain
n = 50;           % Number of points on domain
% BUILD DISCRETE GRID -------------------------------------------------- %
x = linspace(a,b,n); % Set of equally spaced points
h = x(2) - x(1);     % Difference between points
[X,Y] = meshgrid(x); % Creates points to be used for phi function    
% ---------------------------------------------------------------------- %
% LEVEL SET FUNCTION --------------------------------------------------- %
P = ((X).^2 + (2*Y).^2).^(1/2) - 4;
% ---------------------------------------------------------------------- %
% DYNAMIC INTERFACE TRACKING ------------------------------------------- %
figure('units','normalized','outerposition',[0 0 1 1]) % Full screen plot
for i = 1:S
    G = -P;
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
    Pxx = diff(P([1 1:end end],:),2)/h^2;
    Pyy = diff(P(:,[1 1:end end])',2)'/h^2;
    Px = (P(3:end,:)-P(1:end-2,:))/(2*h); 
    Px = Px([1 1:end end],:);
    Py = (P(:,3:end)-P(:,1:end-2))/(2*h); 
    Py = Py(:,[1 1:end end]);
    Pxy = (Px(:,3:end)-Px(:,1:end-2))/(2*h); 
    Pxy = Pxy(:,[1 1:end end]);
    F = (Pxx.*Py.^2-2*Px.*Py.*Pxy+Pyy.*Px.^2)./(Px.^2+Py.^2).^1.5;
    F = -F/100;
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