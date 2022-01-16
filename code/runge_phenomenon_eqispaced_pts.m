% RUNGE PHENOMENON: equispaced collocation pts (Figure 3.6 (a))

myColor = [1 0 1; 1 0 0; 0 0 1; 0 1 0; 0 1 1]; % Color palette

f = @(x) 1./(1+25*x.^2);            % Runge function

x = linspace(-1,1,500);             % Equispaces collocation pts
y=f(x);                             % Eval Runge funtion at collocation pts
plot(x,y,'k','linewidth',2);        % Plot Runge function
hold on;

for N = 3:2:11                      % N = 3,5,7,9,11
    xdata = linspace(-1,1,N+1);     % N+1 equispaced collocation pts
    ydata= f(xdata);                % Eval data
    p = polyfit(xdata,ydata,N);     % Coefficients for polynomial
    y_fit = polyval(p,x);           % Eval
    plot(x,y_fit,'-','Linewidth',1,'Color', myColor((N-1)/2,:));%Plot aprox
end

xlabel('x'), ylabel('f(x)');
legend('Runge function','N=3','N=5','N=7','N=9','N=11');
title('Runge pheomenon. Equispaced collocation pts.');
