% SOLVE LANE-EMDEN EQUATION WITH CHEBFUN PACKAGE


%Color palette
myColor = [1 0 1; 0 0 1; 0 1 0; 0.9290, 0.6940, 0.1250; 1 0 0; 0 0 0];



for n = 0:5                                      % Polytropic index
    N = chebop(@(x,u) x.*diff(u,2)+2.*diff(u)+x.*u.^n,[0,20]);%L-E operator
    N.lbc = @(u) [u-1; diff(u)];                 % Left boundary conditions
    u = solvebvp(N,0);                           % Solve BVP
    plot(u,'Linewidth',1.5,'Color',myColor(n+1,:)), hold on; % Plot 
end


axis([0 20 -0.25 1]);
yline(0);  
xline(0); 
xlabel('\xi'), ylabel('\theta_{n}'); 
title('Lane-Emden solution with Chebfun package for n=0,1,2,3,4,5');
legend('n=0','n=1','n=2','n=3','n=4','n=5');