% EXAMPLE 4: solve BVP u_xx = exp(u), u(-1)=u(1)=0

% Based on Trefethen, L. N. (2000). Spectral methods in MATLAB. Society 
% for industrial and applied mathematics. 


N = 16;                            % No. of interpolation points
[D,x] = cheb(N);                   % Cheb. differentiation matrix & grid
D2 = cheb2(N);                     % Cheb. 2nd order differentiation matrix
D2 = D2(2:N,2:N);                  % Boundary conditions
u = zeros(N-1,1);
change = 1; iter = 0;              % Initialize
while change > 1e-15               % Fixed-point iteration
    unew = D2\exp(u);              % Solve system
    change = norm(unew-u, inf);    % Calculate error
    u = unew; iter = iter + 1;     % Update
end

u = [0;u;0];
clf, subplot('position', [.1 .4 .8 .5]);
plot(x,u,'.','markersize',16); 
xlabel('x');
ylabel('u');
xx = -1:.01:1;
uu = polyval(polyfit(x,u,N),xx);   % Interpolate grid data
line(xx,uu);
grid on;
title(sprintf('no steps = %d    u(0)=%18.14f',iter,u(N/2+1)))
