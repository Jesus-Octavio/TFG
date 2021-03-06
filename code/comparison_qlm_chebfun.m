% comparison QLM-Chebfun
% Boyd, J. P. (2011). Chebyshev Spectral Methods and the Lane-Emden problem.
% Numerical Mathematics: Theory, Methods and Applications, 4(2), 142-157.

format longg                       % Variable precision

N = 90;                            % Number of collocation pts
n = 1;                             % Polytropic index
iters = 10;                        % Number Newton iterations
 
[xi, Y, XCheb, a] = qlm(N,n,iters);
a;
R = xi.*XCheb;                     % x-data
plot(R, Y, 'b-*')                  % Plot function
hold on
plot(R(1), Y(1), 'r-o')            % Mark last point. eigenvalue & 1st root
hold on

N = chebop(@(x,u) x.*diff(u,2)+2.*diff(u)+x.*u.^n,[0,5]);%L-E operator
N.lbc = @(u) [u-1; diff(u)];                 % Left boundary conditions
u = solvebvp(N,0);                           % Solve BVP
plot(u,'Linewidth',1.5), hold on; % Plot 

text(R(1)-0.2,Y(1)+0.1,['\xi=',num2str(R(1)),'...'])
 
disp('First root & eigenvalue xi =  ');
disp(xi);

% fprintf('EL primer cero de la funci?n es:  %.2f.\n', xi); 
% But not format longg

met = 'QSLM';
s=sprintf('Solve Lane-Emden with met=%s, Nodes N=%g, Polytropic index n=%g, Iterations=%g',met,N,n,iters);
title({' ','QSLM. Boyd (SLM)',s,});




