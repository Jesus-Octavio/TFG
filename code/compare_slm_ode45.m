% Compare solutions to the Lane-Emden equation with SLM and ode45

N = 40;                                              % No. of Chebyshev pts
n = 1;                                               % Polytropic index
iters = 8;                                           % No. of iterations

[Yz, ALPHAS, T] = slm(40,1,8);

fprintf('First zero (alpha) is %s', sum(ALPHAS))
figure(1)
Y = sum(Yz,2);
x_values = sum(ALPHAS).*T;

hold on
plot(x_values, Y, 'r-*');

met = 'SLM';
s=sprintf('Solve Lane-Emden with ,met=%s, Chebyshev pts N=%g,polytropic index n=%g,Iterations=%g',met,N,n,iters);
title({' ','Successive linearization method (SLM)',s,})


% Solve Lane-Emden with ode45 for polytropic index n = 0,1,2,3,4,5
options=odeset('AbsTol',0.0000000001,'RelTol',0.00000001,'InitialStep',...
    0.001,'Refine',100);
y0=[1 0];
L=[0.00001 3.2];
n=1;
[T,Y]=ode45(@lane_emden,L,y0,options,n);

plot(T,Y(:,1),'b-*')
grid on
hold off

xL = xlim;
line(xL, [0 0]); 
legend('slm','ode45')
