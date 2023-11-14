clc; clear; close all

N = 10^3;

a = 1.96;


X = rand(N, 1);
Y = rand(N, 1);
W = rand(N, 1);
T = rand(N, 1);

% f(x)
J = 1/N * sum(exp(X.*Y.*W.*T));
Jtheo = 1 + 0.5^4 + 0.5*(1/3)^4 + 1/6*(1/4)^4 + 1/24*(1/5)^4; 
s1 = sqrt(1/(N-1) * (sum(exp(X.*Y.*W.*T).^2) - N*J^2));

% g(x) partie 1
Jtheo2 = 0.5^4;
J2 = 1/N * sum(exp(X.*Y.*W.*T) - X.*Y.*W.*T) + Jtheo2; 
s2 = sqrt( 1/(N-1) * ( sum((exp(X.*Y.*W.*T) - X.*Y.*W.*T + Jtheo2).^2 ) - N*J2^2 ) );

% g(x) partie 2
Jtheo3 = 0.5^4 + 0.5*(1/3)^4; 
J3 = 1/(N) * sum(exp(X.*Y.*W.*T) - X.*Y.*W.*T - 0.5*(X.*Y.*W.*T).^2) + Jtheo3;
s3 = sqrt(1/(N-1) * (sum((exp(X.*Y.*W.*T) - X.*Y.*W.*T - 0.5*(X.*Y.*W.*T).^2 + Jtheo3).^2) - N*J3^2));


intervalle1 = [-a*s1/sqrt(N) + J, a*s1/sqrt(N) + J];
intervalle2 = [-a*s2/sqrt(N) + J2, a*s2/sqrt(N) + J2];
intervalle3 = [-a*s3/sqrt(N) + J, a*s3/sqrt(N) + J];

if intervalle1(1) < J && J < intervalle1(2)
    calculOk = 1;
else
    calculOK = 0;
end