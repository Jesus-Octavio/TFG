% ALIASING PHENOMENON (Theorem 10, Figure 3.3)

% Plot Chebyshev polynomials T_{n} for n = 1,5,7,11,13 and common values

myColor = [1 0 1; 0 0 1; 0 1 0; 1 1 0; 1 0 0; 0 0 0]; % Color palette

N=3;                                           % No. of Chebyshev pts.

i = 0:3;                                       % Compute Chebyshev pts.
z = cos(i*pi/N);                               

for n = [1 5 7 11 13]
    p = poly(chebpoly(n));                     % Chebyshev polynomial
    points = polyval(p ,z);                    % Evaluate at Chebyshev pts.
    plot(chebpoly(n),'Linewidth',1); hold on   % Plot Chebyshev polynomial
end

plot(z, points, 'b-o', 'Markersize', 6, 'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 0 0])              % Plot values at Chebyshev pts.

title('Aliasing. Chebyshev polynomials T_{n}(x), n = 1,5,7,11,13');
xlabel('x');
ylabel('T_{n}(x)');
legend('n=1','n=5','n=7','n=11','n=13');
