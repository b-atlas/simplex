clear all
clc

disp('La méthode du Simplexe ordinaire: \n')
A=input('Donner la matrice A:\n')
b=input('Donner le vecteur b:\n')
c=input('Donner le vecteur c:\n')
[M,x,z] = simOrd(A,b,c);