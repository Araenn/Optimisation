clc; clear; close all

nb_estimation = 10;
erreur = zeros(nb_estimation, 1);
for i = 1:nb_estimation
    erreur(i) = tp1;
end

figure(3)
plot((1:nb_estimation)',erreur)
grid on
title("Erreurs selon les estimations")
xlabel("n estimations")

function erreur = tp1
    t = (0:0.01:1)';
    data = zeros(length(t), 2);
    
    % valeurs à trouver
    x1 = -4;
    x2 = -1;
    x3 = 4;
    x4 = -5;

    bruit = 0.005;
    xv = [x1; x2; x3; x4];
    data(:,1) = t';
    % données
    N = length(xv);
    ym = 0;
    for i = 1:N/2
        ym = ym + xv(N/2 + i) * exp(xv(i) * t);
    end
    ym = ym + bruit * randn(length(t), 1);
    % hm = ym - xv(3) * exp(xv(1) * t) + xv(4) * exp(xv(2) * t)+ bruit * randn(length(t), 1);
    data(:,2) = ym;

    % initialisation
    x0 = [-5;-10;-4;2];
    x0 = randn(N, 1);
    niter = 1000;
    critere = 10^-12;

    % premier jacobien
    J1 = zeros(length(t), N/2);
    J2 = J1;
    for i = 1:N/2
        J1(:, i) = -x0(N/2+i)*t.*exp(x0(i)*t);
        J2(:, i) = -exp(x0(i)*t);
    end
    Ju0 = [J1, J2];
    %premier hessien
    Hu0 = Ju0' * Ju0;
    %premier mu
    Hu0 = inv(Hu0);

    [x_estim, f] = BFGS(niter, critere, data, x0, Hu0);
    if isnan(x_estim)
        x_estim = zeros(N, 1);
    end
    x_estim
    figure(1)
    semilogy((0:niter-1)', f)
    hold on
    grid on
    title("Convergence de la fonction de coût")
    xlabel("Nombre d'itérations")
    ylabel("Fonction de coût")

    y_est = 0;
    for i = 1:N/2
        y_est = y_est + x_estim(N/2 + i) .* exp(x_estim(i) * t);
    end
    
    erreur = mean(ym.^2-y_est.^2);
    figure(2)
    plot(t, ym, 'b')
    hold on
    plot(t, y_est, 'o')
    grid on
    legend("Données", "Estimation")
    title("Analyse de l'estimation")
end

function [x_estim, f] = BFGS(niter, critere, data, x0, Hu0)
    
    f = zeros(niter, 1);
    f(1) = 10^10;
    erreur = 10;
    i = 2;
    
    % Matrice BFGS initiale
    B = Hu0;
    
    while erreur > critere
        if i > niter
            break;
        end
        [f_courant, g_courant, J] = fcout(data, x0);
        
        
        % Direction de recherche
        d = -B \ g_courant;
        
        % Calcul du pas (alpha) avec une recherche linéaire
        alpha = linesearch(data, x0, d, 0.0001);

        % Mise à jour de x
        x_next = x0 + alpha * d;
        
        [f_next, g_next, J_next] = fcout(data, x_next);
        
        y = g_next - g_courant;
        s = alpha * d;
        z = B*s;

        % Mise à jour de la matrice BFGS
        B = B + (y*y')/(y'*s) - (z*z')/(s'*z);
        
        f(i) = f_next;
        erreur = abs((f(i) - f(i-1)))/f(i-1);
        
        
        
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
    N = length(x);
    % calcul hm
    h = data(:, 2);
    for i = 1:N/2
        h = h - x(N/2 + i) .* exp(x(i) * t);
    end
    % fonction cout
    f_courant = (1/2) * norm(h)^2;
    J1 = zeros(length(t), N/2);
    J2 = J1;
    for i = 1:N/2
        J1(:, i) = -x(N/2+i)*t.*exp(x(i)*t);
        J2(:, i) = -exp(x(i)*t);
    end
    J = [J1, J2];
    G = J' * h;
end