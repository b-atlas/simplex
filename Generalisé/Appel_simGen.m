clear all
clc

disp('La m�thode du Simplexe generalisee: \n')
A=input('Donner la matrice A:\n')
b=input('Donner le vecteur b:\n')
chr=input('Donner le string chr:\n')
    %correspond a saisir la nature in�galit�s et �galit�s dans l'ordre
    %pour entrer <= il faut taper <
    %pour entrer >= il faut taper >
    %pour entrer = il faut taper =
c=input('Donner le vecteur c:\n')
[nv_M,x,z] = simGen(A,chr,b,c);