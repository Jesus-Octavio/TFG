% SOLVE LANE-EMDEN EQUATION WITH ode45

options=odeset('AbsTol',0.0000000001,'RelTol',0.00000001,'InitialStep',...
    0.001,'Refine',100);

myColor = [1 0 1; 0 0 1; 0 1 0; 1 1 0; 1 0 0; 0 0 0]; % Color palette

for n = 0:5                % Polytropic index
    y0=[1 0];              % Initial condictions
    L=[0.00001 5];         % epsilon = 0.00001 (Avoid singularity at oigin)
    [T,Y]=ode45(@lane_emden,L,y0,options,n);
    plot(T,Y(:,1),'Linewidth',1.5,'Color',myColor(n+1,:)), hold on;
end

axis([0 5 -1 1]);
yline(0);  
xline(0); 
xlabel('\xi'), ylabel('\theta'); 
title('Lane-Emden equation with ode45 for n=0,1,2,3,4,5');
legend('n=0','n=1','n=2','n=3','n=4','n=5');
