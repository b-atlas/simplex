clear all 
clc

disp('            La méthode du Simplexe ordinaire: ')
A=input('Donner la matrice A:\n')
b=input('Donner le vecteur b:\n')
c=input('Donner le vecteur c:\n')
format shortG % pour convertir les valeurs du tableau de 1.0e+04 en normal 
[Tab,x,z]=simOrd(A,b,c)
