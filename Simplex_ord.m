clear all;
close all;
clc;

function z = simOrd(A, b, z)

[n1,m1] = size(A);
M = [A, eye(n1), b];
[n2,m2] = size(M);
z = [z, zeros(1, m2 - length(z) + 1)];
M = [M; z];
disp(M)
end