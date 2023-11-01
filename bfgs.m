clc; clear; close all

tp1

function tp1
    t = (0:0.01:1)';
    data = zeros(length(t), 2);
    
    % valeurs à trouver
    x1 = -4;
    x2 = -1;
    x3 = 4;
    x4 = -5;

    bruit = 0;
    xv = [x1; x2; x3; x4];
    data(:,1) = t';

    % données
    ym = xv(3) * exp(xv(1) * t) + xv(4) * exp(xv(2) * t)+ bruit * randn(length(t), 1);
    data(:,2) = ym;

    % initialisation
    x0 = [-5;-10;-4;2];
    x0 = randn(4, 1)
    niter = 1000;
    critere = 10^-12;

    % premier jacobien
    Ju0 = [-x0(3)*t.*exp(x0(1) * t),  - x0(4)*t .* exp(x0(2)*t), -exp(x0(1)*t), - exp(x0(2)*t)];
    %premier hessien
    Hu0 = Ju0' * Ju0;
    Hu0 = inv(Hu0);

    [x_estim, f] = BFGS(niter, critere, data, x0, Hu0);
    x_estim

    figure(1)
    semilogy((0:niter-1)', f)
    grid()
    title("Convergence de la fonction de coût")
    xlabel("Nombre d'itérations")
    ylabel("Fonction de coût")
end

function [x_estim, f] = BFGS(niter, critere, data, x0, Hu0)
    
    f = zeros(niter, 1);
    f(1) = 10^10;
    erreur = 10;
    i = 2;
    
    % Matrice BFGS initiale
    H = Hu0;
    
    while erreur > critere
        [f_courant, g_courant, J] = fcout(data, x0);
        
        
        % Direction de recherche
        d = -H \ g_courant;
        
        % Calcul du pas (alpha) avec une recherche linéaire
        alpha = linesearch(data, x0, d, 0.0001);

        % Mise à jour de x
        x_next = x0 + alpha * d;
        
        [f_next, g_next, J_next] = fcout(data, x_next);
        
        y = g_next - g_courant;
        s = alpha * d;
        z = H*s;

        % Mise à jour de la matrice BFGS
        H = H + (y*y')/(y'*s) - (z*z')/(s'*z);
        
        f(i) = f_next;
        erreur = abs((f(i) - f(i-1)))/f(i-1);
        
        if i > niter
            break;
        end
        
        x0 = x_next; % Mise à jour de x0
        i = i + 1;
    end
    
    x_estim = x0;
end


function alpha = linesearch(data, x0, d, tol)
    f = @(alpha) fcout(data, x0 + alpha * d);  % fonction coût en fonction de alpha

    % minimum de la fonction de coût en fonction de alpha
    alpha = fminbnd(f, 0, 1, optimset('TolX', tol));
end


function [f_courant, G, J] = fcout(data, x)
    t = data(:, 1);

    % calcul hm
    h = data(:, 2) - x(3) .* exp(x(1) * t) - x(4) .* exp(x(2)*t);

    % fonction cout
    f_courant = (1/2) * norm(h)^2;
    J = [-x(3)*t.*exp(x(1) * t),  - x(4)*t .* exp(x(2)*t), -exp(x(1)*t), - exp(x(2)*t)];
    G = J' * h;
end