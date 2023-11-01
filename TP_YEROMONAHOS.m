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
    % hm = ym - xv(3) * exp(xv(1) * t) + xv(4) * exp(xv(2) * t)+ bruit * randn(length(t), 1);
    data(:,2) = ym;

    % initialisation
    x0 = [-5;-10;-4;2];
    niter = 1000;
    critere = 10^-12;
    tau = 1^-6;

    % premier jacobien
    Ju0 = [-x0(3)*t.*exp(x0(1) * t),  - x0(4)*t .* exp(x0(2)*t), -exp(x0(1)*t), - exp(x0(2)*t)];
    %premier hessien
    Hu0 = Ju0' * Ju0;
    %premier mu
    mu0 = tau * max(diag(Hu0));
    nu = 2;

    [x_estim, f] = LM(niter, critere, mu0, data, x0, nu);
    x_estim
    figure(1)
    semilogy((0:niter-1)', f)
    grid()
    title("Convergence de la fonction de coût")
    xlabel("Nombre d'itérations")
    ylabel("Fonction de coût")
end

function [x_estim, f] = LM(niter, critere, mu0, data, x0, nu)
    
    f = zeros(niter, 1);
    f(1) = 10^10;
    erreur = 10;
    i = 2;
    while erreur > critere
        [x_estim, f_courant, mu_courant, nu_courant] = dir_LM(x0, data, mu0, nu);

        f(i) = f_courant;
        erreur = abs((f(i) - f(i-1)))/f(i-1);

        if i > niter
            break;
        end

        mu0 = mu_courant;
        nu = nu_courant;
        x0 = x_estim;
        
        
        i = i + 1;
    end
    x_estim = x0;
end

function [x_estim, f_courant, mu_courant, nu_courant] = dir_LM(x0, data, mu0, nu)
   x_estim = zeros(4, 1);
   % premiere fonction cout
   [f, G, J] = fcout(data, x0);
   
   %direction descente
   dLM = -inv(J' * J + mu0*eye(4)) * G;
   
   % deuxieme fonction cout
   xa = x0 + dLM;
   [f_courant, G2, J2] = fcout(data, xa);

   % calcul gamma
   u = f'*f;
   u_dLM =  0.5 * u + dLM'*G + 0.5 * dLM' * (J') * J * dLM;

   gamma = (f + f_courant) / (u - u_dLM);
   
   % heuristique
   if gamma > 0
       x_estim = x0 + dLM;
       mu_courant = mu0 * max(1/3, 1 - (2 * gamma - 1)^3);
       nu_courant = 2;
   else
       mu_courant = mu0 * nu;
       nu_courant = 2 * nu;
   end
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