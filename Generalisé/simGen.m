function [nv_M,x,z] = simGen(A,chr,b,c)

%Début de phase 1

%construction du tableau initial de phase 1 
k=1;
nb_iterations = 20;
[n1,m1] = size(A);
identite = eye(n1);

nbr_e = 0;

M = A;

% on ajoute ici les colonnes des variables e1, e2 ..
for i = 1:n1
   if(chr(i) == '>')
       M = [M -identite(:,i)];
       nbr_e = nbr_e+1;
   elseif(chr(i) == '<')
       M = [M identite(:,i)];
       nbr_e = nbr_e+1;
   end
end

side = []; % un vecteur qui va contenir les variables de bases
for i = 1:n1 %initialiser les variables de base
    side(i) = m1+nbr_e+i;
end

M = [M eye(n1) b]; %on ajoute les colonnes des variables y en plus du b

[n2,m2] = size(M);
z_pr = zeros(1,m2); %le vecteur z'

for j = 1:m2 %on affecte à chaque element de z' sa valeur correspondant
   if (j<=m1+nbr_e | j==m2)
       for i = 1:n2
           z_pr(j) = z_pr(j) + M(i,j);
       end
   end
end

M = [M; z_pr]; %on ajoute le vecteur ligne z' à la fin de la matrice
[n2,m2] = size(M);

%fin de construction du tableau initial de phase 1
disp('Le tableau initial de la phase 1 est:')
disp(M)

%début du simplex ordinaire de phase 1
for iteration = 1:nb_iterations
    tab_fin=M(end,1:end-1)>0; %retourne 0 si tous les elements de tab_fin <= 0
    if(tab_fin==0) %condition d'arret
        break;
    end;
    
    %determiner le pivot
    [cpp,col_piv] = max(M(end,1:end-1));
    rap = M(1:end-1,end)./M(1:end-1,col_piv); %determiner les rapports
    
    indices = rap<0; %on donne aux rapports négatifs 
    rap(indices) = Inf; %la valeur infini
    
    [ppr,lin_piv] = min(rap);
    pivot = M(lin_piv,col_piv); %le pivot est determiné
    
    side(lin_piv) = col_piv; %on remplace la variable sortante par la variable entrante dans les variables de base
    
    for i = 1:n2 %on effectue les operations necessaires sur les lignes de la matrice (sauf la ligne de pivot)
        if(i ~= lin_piv)
            x = M(i,col_piv)/pivot;
            M(i,:) = M(i,:)-x*M(lin_piv,:);
        end
    end
    
    M(lin_piv,:)= M(lin_piv,:)/pivot; %on effectue l'opération necessaire sur la ligne de pivot
    
end
%fin du simplex ordinaire de phase 1

if(abs(M(end,end))<1e-14) %juste pour fixer un probleme de précision au programme
   M(end,end)=0; 
end

disp('Le tableau final de la phase 1 est:')
disp(M);

%début de phase 2

if(tab_fin == 0 & M(end,end)==0) %tester si on peut passer de phase 1 à phase 2
    
    for i=1:n1 %encore on teste si on peut passe de phase 1 à phase 2 (sinon on arrete la fonction)
        if (side(i)>m1+nbr_e & side(i)<m2 & M(i,end)~=0)
            disp('Le programme linéaire ne posséde aucune solution de base admissible1!')
            %mettre ces affectations est necessaire juste pour des questions
            %d'implémentaion, les outputs doivent forcement avoir des valeurs
            nv_M = M;
            z = nv_M(end,end);
            x=0;
            return
        end
    end
    
    %construction du tableau initial de phase 2
    nv_M = [M(1:end-1,1:m1+nbr_e)]; 
    for i = 1:n1
        if(side(i) > m1+nbr_e & side(i)<m2)
            nv_M = [nv_M M(1:end-1,side(i))];
        end
    end
    nv_M = [nv_M M(1:end-1,end)];
    [n3,m3] = size(nv_M);
    
    %début construction de z
    z=zeros(1,m3);
    deBase=zeros(1,m3-1); %vecteur qui contient 1 pour l'indice de vecteur de base, sinon contient 0
    for i = 1:m3-1 %on remplit ce vcteur
        for j=1:n1
            if(side(j)==i)
                deBase(i)=1; 
            end
        end
    end
    deBase(m3)=0;
    indice1=1;
    
    for i=1:m3
        if(deBase(i)==1) %si la vecteur de base on met 0 au case correspondant
           z(i)=0; 
        else
            for j=1:i %on parcour les colonnes de 1 à i
                if(deBase(j)==1) %chechant les variables de bases
                    for h=1:n1 %on trouve l'indice de valeur 1 dans la colonne
                       if j==side(h)
                          indice1=h;
                       end
                    end
                    
                    if j<=numel(c)
                        z(i) = z(i)+c(j)*nv_M(indice1,i); %et on fait l'addition entre z(i) et la multiplication de coefficient corespondant dans z
                           % avec l'element de la matrice de meme ligne
                           % que 1 et le meme col que l'elelement de z
                    end
                end
                
            end
        end
    end
    
    %fin de construction de z
    
    nv_M = [nv_M; z];
    [n3,m3] = size(nv_M);
    
    %fin de construction du tableau initial de phase 2
    disp('Le tableau initial de la phase 2 est:')
    disp(nv_M)
    
    %simplex ordinaire de phase 2
    for iteration = 1:nb_iterations
        tab_fin=nv_M(end,1:end-1)>0;
        if(tab_fin==0)
            break;
        end;

        [cpp,col_piv] = max(nv_M(end,1:end-1));
        rap = nv_M(1:end-1,end)./nv_M(1:end-1,col_piv);

        indices = rap<0;
        rap(indices) = Inf;

        [ppr,lin_piv] = min(rap);
        pivot = nv_M(lin_piv,col_piv);

        for i = 1:n3
            if(i ~= lin_piv)
                x = nv_M(i,col_piv)/pivot;
                nv_M(i,:) = nv_M(i,:)-x*nv_M(lin_piv,:);
            end
        end

        nv_M(lin_piv,:)= nv_M(lin_piv,:)/pivot;
    end
    
    %fin de phase 2
    
    disp('Le tableau final de la phase 2 est:')
    disp(nv_M)
    
    %stockage des résultats
    z=-nv_M(end,end);
    
    identite = eye(n3);
    for i = 1:n3 %determiner les solutions x
        x(i)=0;
        for j = 1:n3
           if(nv_M(:,i) == identite(:,j))
              x(i) = nv_M(j,end); 
           end
        end
    end
  
else %si on peut pas passer de phase 1 à 2 on affiche ce message
    disp('Le programme linéaire ne posséde aucune solution de base admissible.')
    %mettre ces affectations est necessaire juste pour des questions
    %d'implémentaion, les outputs doivent forcement avoir des valeurs
    nv_M = M;
    z = nv_M(end,end); 
    x=0;
end

end

