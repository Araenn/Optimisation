clc; clear; close all

n = [10^1, 10^2, 10^3, 10^4, 10^5, 10^6]';

a = 1.96;

% tirer un point au hasard, calculer toutes distances et determiner plus
% courte, puis determiner nouveau point tel que segment entre les deux
% (plus qu'un point apres), refaire bcp de fois
% doit ressembler à une mosaique si marche bien

for i = 1:length(n)
    N = n(i);
    XMC = rand(N, 1);
    YMC = rand(N, 1);
    
end

figure(1)
affichage = 100;
scatter(XMC(1:affichage), YMC(1:affichage))
grid()
