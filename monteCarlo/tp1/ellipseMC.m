clc; clear; close all

N = 10^3;

A = 1/2;
B = 1/2;
a = 1.96;

X = ((rand(N, 1) - A) * 2).^2;
Y = ((rand(N, 1) - B) * 2).^2;

Z = (X + Y) < 1;

n1 = sum(Z);

p = n1/N;
s = sqrt(p * (1-p));
S = 4*A*B*p;
Stheo = A*B*pi;

intervalle = [-a*s*A*B/sqrt(N) + S, a*s*A*B/sqrt(N) + S];

if intervalle(1) < p && p < intervalle(2)
    calculOk = 1;
else
    calculOK = 0;
end