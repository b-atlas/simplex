function [Tab,x,z]=simOrd(A,b,c)

nb_iter=5; %le nombre d'itération
n=size(A,1); %nombre de ligne de la matrice A
m=size(A,2); %nombre de colonne de la matrice A
Tab=[A eye(n) b;c zeros(1,n+1)]; %Le tableau initial
x=0; %le point optimal
z=0; %la valeur de la fonction objectif
a=0;

for iter=1:nb_iter
    Tab_fin=Tab(end,1:end-1)>0; %retourne 1 s'il existe des valeurs positive dans la ligne de z 
    if(Tab_fin==0)% si la condition est fausse alors on atteint la valeur optimale de z
        break; %on sort de la boucle
    end
    
    [cpp,col_pivot]=max(Tab(end,1:end-1)) %cpp: coef plus positive (la valeur),col_pivot: l'indice de col de cpp 
    Rap=Tab(1:end-1,end)./Tab(1:end-1,col_pivot) %le rapport: on divise les valeurs de b sur les valeur d col_pivot 
    ind1=Rap<=0 %stock les rapport negatif
    Rap(ind1)=Inf; %on affecte l'infinie au rapports negatif
    [mn,ligne_pivot]=min(Rap) %min: minimum des rapports, ligne_pivot: l'indice de la ligne de pivaot
    pivot=Tab(ligne_pivot,col_pivot); %le pivot est l'intercection de la ligne de pivot et de la colonne de pivot
    for i=1:n+1 % une boucle de la premier ligne de tableau jusqu'à la dernier ligne (z)
        if(Tab(i,col_pivot)~=Tab(ligne_pivot,col_pivot)) % si la ligne n'est pas la ligne de pivot
           coef=Tab(i,col_pivot); 
           Tab(i,1:end)=Tab(i,1:end)-(coef/pivot)*Tab(ligne_pivot,1:end) % ex: L1-coef/pivot*L2
        end
    end
    Tab(ligne_pivot,1:end)=Tab(ligne_pivot,1:end)./pivot % on divise la ligne de pivot sur le pivot
end

%------------ la valeur de la fonction objectif ---------------

z=-Tab(end,end);

%------------ les coords de x,le point optimale ---------------

%for i=1:m %une boucle sur les colonnes de A (le nombre de colonnes est le nbr de coords de x) 
%    x(i)=0 %si la variable n'est pas une variable de base alors sa coord = 0
 %   for j=1:n % boucle sur les lignes
 %       if(Tab(j,i)==1) %l'indice de la valeur de a est l'intersection de la premiere colonne et la ligne dont la valeur vaut 1
  %          x(i)=Tab(j,end)
   %     end
    %end
%end

%-------------------------------------------------------------
id=eye(n+1);

for i=1:m 
    x(i,1)=0;
    for j=1:n+1
        if Tab(:,i)==id(:,j);
            x(i,1)=Tab(j,end)
        end
    end
end
