clc; clear; close all

n = [10^1, 10^2, 10^3, 10^4, 10^5, 10^6]';

a = 1.96;

% ici on calcule f(x)

for i = 1:length(n)
    N = n(i);
    XMC = rand(N, 1);
    YMC = rand(N, 1);
    
    IMC(i) = 1/N * sum(exp(XMC+YMC));
    
    p = (1:N)';
    XQMC = p * sqrt(2) - fix(p * sqrt(2));
    YQMC = p * sqrt(3) - fix(p * sqrt(3));
    
    IQMC(i) = 1/N * sum(exp(XQMC+YQMC));
end

figure(1)
affichage = 100;
scatter(XMC(1:affichage), YMC(1:affichage))
hold on
scatter(XQMC(1:affichage), YQMC(1:affichage))
grid()
legend("Monte-Carlo", "Quasi MC")