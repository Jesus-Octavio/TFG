close all
N = 100; % number of chebyshev nodes
iters = 5; % iterations

myColor = [1 0 1; 0 0 1; 0 1 0; 1 1 0; 1 0 0; 0 0 0];
for n = 0:5 %polytropic index
    [Yz, xis, T] = slm(N,n,iters);
    Y = sum(Yz,2);
    x_values = sum(xis).*T;
    plot(x_values,Y,'Linewidth',1.5,'Color',myColor(n+1,:)), hold on;
    grid on    
    axis([0 5 -1 1]);
    title('Solution of the Lane-Emden equation with slm for n=0,1,2,3,4,5');
    xlabel('t'), ylabel('x')
end


% 
xL = xlim;
yL = ylim;
line([0 0], yL);  %y-axis
line(xL, [0 0]);  %x-axis   
legend('n=0','n=1','n=2','n=3','n=4','n=5');

% print -dpng  poly-slm-todo.png 


