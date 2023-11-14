clc; clear; close all

N = 10^3;

A = 1/2;
B = 1/2;
r = A^2;
a = 1.96;

X = (rand(N, 1) - A).^2;
Y = (rand(N, 1) - B).^2;

Z = (X + Y) < r;

n1 = sum(Z);

p = n1/N;
s = sqrt(p * (1-p));
S = A*B*p;
Stheo = pi*r^2;

intervalle = [-a*s/sqrt(N) + p, a*s/sqrt(N) + p];

if intervalle(1) < p && p < intervalle(2)
    calculOk = 1;
else
    calculOK = 0;
end