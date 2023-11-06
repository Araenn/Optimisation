clc; clear; close all

nb_estimation = 1;
erreur = zeros(nb_estimation, 1);
for i = 1:nb_estimation
    erreur(i) = tp1;
end

% figure(3)
% stem((1:nb_estimation)',abs(erreur))
% grid on
% title("Erreurs selon les estimations")
% xlabel("n estimations")

function erreur = tp1
    t = (0:0.01:10)';
    data = zeros(length(t), 2);
    
    % valeurs à trouver

    bruit = 0;
    xv = [-2; -1.5; -1; -0.5; -0.2; 10; 5; 2; -5; -10];
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
    x0 = randn(N, 1);
    x0 = [0.8404   -0.8880    0.1001   -0.5445    0.3035   -0.6003    0.4900    0.7394    1.7119   -0.1941]';
    %x0 = [-0.3899; 1.3185; 1.8374; 0.8807; -0.5136; -1.4283; -0.3721; -1.0776; -0.3248; -0.3665]; % non
    %x0 = [1.3054; 0.3658; -1.5758; 1.0079; 0.4747; -1.5317; -0.0272; -0.4168; -1.7485; -0.4404]; %non
    niter = 1000;
    critere = 10^-12;
    tau = 1^-6;

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
    mu0 = tau * max((Hu0(:)));
    nu = 2;

    [x_estim, f] = LM(niter, critere, mu0, data, x0, nu);
    
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

function [x_estim, f] = LM(niter, critere, mu0, data, x0, nu)
    
    f = zeros(niter, 1);
    f(1) = 10^10;
    erreur = 10; % plus d'iterations
    i = 2;
    while i < niter
        if i > niter
            break;
        end
        [x_estim, f_courant, mu_courant, nu_courant] = dir_LM(x0, data, mu0, nu);

        f(i) = f_courant;
        erreur = abs((f(i) - f(i-1)))/f(i-1);

        mu0 = mu_courant;
        nu = nu_courant;
        x0 = x_estim;
        
        
        i = i + 1;
    end
    x_estim = x0;
end

function [x_estim, f_courant, mu_courant, nu_courant] = dir_LM(x0, data, mu0, nu)
   N = length(x0);
   % premiere fonction cout
   [f, G, J] = fcout(data, x0);
   
   %direction descente
   dLM = -inv(J' * J + mu0*eye(N)) * G;
   
   % deuxieme fonction cout
   xa = x0 + dLM;
   [f_courant, G2, J2] = fcout(data, xa);

   % calcul gamma

   u_dLM =  f + dLM'*G + 0.5 * dLM' * (J') * J * dLM;

   gamma = (f - f_courant) / (f - u_dLM);
   
   % heuristique
   if gamma > 0
       x_estim = x0 + dLM;
       mu_courant = mu0 * max(1/3, 1 - (2 * gamma - 1)^3);
       nu_courant = 2;
   else
       mu_courant = mu0 * nu;
       nu_courant = 2 * nu;
       f_courant = f;
       x_estim = x0;
   end
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