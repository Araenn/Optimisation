clc; clear; close all

n = 100;
N = 100000;

% tirer un point au hasard, calculer toutes distances et determiner plus
% courte, puis determiner nouveau point tel que segment entre les deux
% refaire bcp de fois
% doit ressembler à une mosaique si marche bien

XMC = rand(n, 1);
YMC = rand(n, 1);
affX = XMC;
affY = YMC;
for i = 1:N

    randX = rand(1);
    randY = rand(1);
    epsilon = 0.01 - (i/N)*(0.01 - 0.0001);

    for j = 1:n
        distance(j) = ((XMC(j)-randX)^2 + (YMC(j)-randY)^2);
    end
    [~, indiceMinDist] = min(distance);
    newX = (1-epsilon)*XMC(indiceMinDist) + epsilon*randX;
    newY = (1-epsilon)*YMC(indiceMinDist) + epsilon*randY;
    
    
    XMC(indiceMinDist) = newX;
    YMC(indiceMinDist) = newY;
end

figure(1)
affichage = n;
scatter(affX(1:affichage), affY(1:affichage))
grid()
hold on
scatter(XMC(1:affichage), YMC(1:affichage))
legend("Anciennes sources", "Nouvelles sources")