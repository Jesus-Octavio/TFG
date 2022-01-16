% CHEBYSHEV POINTS: GEOMETRIC INTERPRETATION

% Chebyshev points ara the projections onto the x-axis of equally spaced
% points on the unit circle. Note they are numberes from right to left.
% (Figure 3.3)

n=16;                                               % No. of Chebyshev pts
tt = linspace(0,pi,n+1);
zz = exp(1i*tt); plot(zz,'.-k', 'Markersize',15);   % Unit circumference
xx = chebpts(n+1);                                  % Chebfun package   
hold on, for j=2:n, plot([xx(n+2-j) zz(j)], 'k'), end
plot(xx,0*xx,'.r', 'Markersize',15)