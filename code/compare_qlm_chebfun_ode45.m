% Interactive program to compare solutions to the Lane-Emden equation
% with different methods: QLM, Chebfun, ode45

% Example:
% User introduces
%   N = 40
%   n = 3
%   iters = 20

format longg 

prompt1 = 'Introduce no. collocation points:  ';
N = input(prompt1);         % No. of collocation pts

prompt2 = 'Introduce polytropic index in [0,5):  ';
n = input(prompt2);          % Polytropic index

prompt3 = 'Introduce number Newton iterations for qlm:  ';
iters = input(prompt3);     % No. of Newton iterations

% Compute QLM method
[xi, Y, XCheb] = qlm(N,n,iters);
R = xi.*XCheb;
hold on
plot(R, Y, 'r-*')

% Compute Chebfun solution
N = chebop(@(x,u) x.*diff(u,2)+2.*diff(u)+x.*u.^n,[0,20]);
N.lbc = @(u) [u-1; diff(u)];
u = solvebvp(N,0);
plot(u, 'g-o'), hold on;

% Compute ode45 solution
options=odeset('AbsTol',0.0000000001,'RelTol',0.00000001,'InitialStep',0.001,'Refine',100);
y0=[1 0];
L=[0.00001 20];
[T,Y]=ode45(@lane_emden,L,y0,options,n);
plot(T,Y(:,1), 'black') %si se marcan los puntos, se tapa todo

xL = xlim;
line(xL, [0 0]);  %x-axis 
hold off
legend('QLM','Chebfun','ode45')