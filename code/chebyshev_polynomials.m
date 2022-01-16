% CHEBYSHEV POLYNOMIALS (Definition 5)

% Plot chebyshev polynomial T_{n} for n = 0,...,4 (Figure 3.2)

myColor = [1 0 1; 0 0 1; 0 1 0; 1 1 0; 1 0 0; 0 0 0]; % Color palette
for N=0:4
    plot(chebpoly(N),'Linewidth',1.5,'Color',myColor(N+1,:)); hold on
end

title('Chebyshev polynomials T_{n}(x) n = 0,...,4');
ylim([-1 1.02]);
xlabel('x');
ylabel('T_{n}(x)');
legend('n=0','n=1','n=2','n=3','n=4');
