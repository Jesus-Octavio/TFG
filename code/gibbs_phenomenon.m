% GIBBS PHENOMENON (Figure 3.5)

% Based on Driscoll, T. A., Hale, N., & Trefethen, L. N. (2014). 
% Chebfun guide.

f = chebfun('sign(x)', 10);                   % Chebfun. 10 collocation pts
subplot (1,2,1), plot(f, '-', 'Linewidth', 2), grid on;
legend('N=10')
xlabel('x');
ylabel('sign(x)');

f = chebfun('sign(x)', 20);                   % Chebfun. 20 collocation pts
subplot(1,2,2), plot(f,'-', 'linewidth', 2), grid on;
legend('N=20');
xlabel('x');
ylabel('sign(x)');

sgtitle('Gibbs phenomenon')
