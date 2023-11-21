clc; clear; close all

n = 9;
N = 50;
M = 10000;

% calcul cellue voronoi = ensemble points plus proches d'une source
% source j remplacee par barycentre cellule voronoi

XMC = rand(n, 1);
YMC = rand(n, 1);
affX = XMC;
affY = YMC;
for i = 1:N

    randX = rand(M, 1);
    randY = rand(M, 1);
    
    for j = 1:n
        for k = 1:M
            distance(j, k) = ((XMC(j)-randX(k))^2 + (YMC(j)-randY(k))^2);
            [~, indiceMinDist(k)] = min(distance(:, k));
            celluleVoronoi(j,k) = (indiceMinDist(k) == j);
        end
        barX(j) = mean(randX(celluleVoronoi(j,:)==1));
        barY(j) = mean(randY(celluleVoronoi(j,:)==1));

        newX(j) = barX(j);
        newY(j) = barY(j);
        
    end
    
    XMC = newX;
    YMC = newY;
    
    
end

figure(1)
affichage = n;
scatter(affX(1:affichage), affY(1:affichage))
grid()
hold on
scatter(XMC(1:affichage), YMC(1:affichage))
legend("Anciennes sources", "Nouvelles sources")