clc; clear; close all

N = 10^3;

a = 1.96;

% ici on calcule f(x)

X = rand(N, 1);
Y = rand(N, 1);
W = rand(N, 1);
T = rand(N, 1);

I = 1/N * sum(exp(X+Y+W+T));
Itheo = (exp(1) - 1)^4;

s = sqrt(1/(N-1) * (sum(exp(X+Y+W+T).^2) - N*I^2));

Stheo = Itheo;

intervalle = [-a*s/sqrt(N) + I, a*s/sqrt(N) + I];

if intervalle(1) < s && s < intervalle(2)
    calculOk = 1;
else
    calculOK = 0;
end