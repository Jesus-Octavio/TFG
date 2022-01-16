% EXAMPLE 3: solve linear BVP u_xx = exp(4x), u(-1)=u(1)=0 

% Based on Trefethen, L. N. (2000). Spectral methods in MATLAB. Society 
% for industrial and applied mathematics. 

N = 16;                            % No. of interpolation points
[D,x] = cheb(N);                   % Cheb. differentiation matrix & grid
D2 = cheb2(N);                     % Cheb. 2nd order differentiation matrix
D2 = D2(2:N,2:N);                  % Boundary conditions
f = exp(4*x(2:N));                 % u_xx = exp(4x)
u = D2\f;                          % Poisson eq. solved here
u = [0;u;0];

clf, subplot('position', [.1 .4 .8 .5]);
plot(x,u,'.','markersize',16);
xlabel('x');
ylabel('u');
xx = -1:.01:1;
uu = polyval(polyfit(x,u,N),xx);   % Interpolate data grid
line(xx,uu);
grid on;
exact=(exp(4*xx)-sinh(4)*xx-cosh(4))/16;   % Analytic solution
title(['max err =' num2str(norm(uu-exact,inf))], 'fontsize', 12)