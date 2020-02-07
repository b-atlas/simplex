function [M,x,z] = simOrd(A)

M = A;
[n2,m2] = size(M);
nb_iterations = 5;

for iteration = 1:nb_iterations
    tab_fin=M(end,1:end-1)>0;
    if(tab_fin == 0)
        break;
    end;

    [cpp,col_piv] = max(M(end,1:end-1));
    rap = M(1:end-1,end)./M(1:end-1,col_piv);
    
    indices = rap<=0;
    rap(indices) = Inf;
    
    [ppr,lin_piv] = min(rap);
    pivot = M(lin_piv,col_piv);

    for i = 1:n2
        if(i ~= lin_piv)
            x = M(i,col_piv)/pivot;
            M(i,:) = M(i,:)-x*M(lin_piv,:);
        end
    end
    
    M(lin_piv,:)= M(lin_piv,:)/pivot;
    
end

z = -M(end, end);

identite = eye(n2);
for i = 1:n2
    x(i)=0;
    for j = 1:n2
       if(M(:,i) == identite(:,j))
          x(i) = M(j,end); 
       end
    end
end

end

