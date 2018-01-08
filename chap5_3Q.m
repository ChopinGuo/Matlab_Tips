clear
close all;

tol = 400*10^(-6);

[np, dp] = pade(tol,6);

delay = tf(np,dp);
delta = tf(np,dp)-1;