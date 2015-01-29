
nconn = size(linelist, 2);

tmp = [linelist.c];

Cent(:,1) = tmp(1:2:end);
Cent(:,2) = tmp(2:2:end);

for i=1:nconn

    %D = [(1:size(Cent,1))' sqrt((Cent(:,1)-linelist(i).c(1)).^2 + (Cent(:,2)-linelist(i).c(2)).^2)];
    
    %D = sortrows(D,2);
    
    nhoodS = Cliques{linelist(i).clq(1)}(Cliques{linelist(i).clq(1)}~=i);
    nhoodE = Cliques{linelist(i).clq(2)}(Cliques{linelist(i).clq(2)}~=i);
    
    if (isempty(nhoodS))
        dS = [];
    else
        dS = [nhoodS' abs([linelist(nhoodS).ang] - linelist(i).ang)'];
        dS = sortrows(dS,2);
    end
    
    if (isempty(nhoodE))
        dE = [];
    else
        dE = [nhoodE' abs([linelist(nhoodE).ang] - linelist(i).ang)'];
        dE = sortrows(dE,2);
    end
    
    if (isempty(dS))
        
        D = [(1:size(Cent,1))' sqrt((Cent(:,1)-linelist(i).c(1)).^2 + (Cent(:,2)-linelist(i).c(2)).^2)];
        D = sortrows(D,2);
        
        nodes = D(:,1);
        nodes = nodes(nodes ~= i);
        
        if (isempty(dE)) 
            G{i} = [i nodes(1) nodes(2)];
        elseif (size(dE,1) == 1)
            if(nodes(1) ~= dE(1,1))
                G{i} = [i nodes(1) dE(1,1)];
            else
                G{i} = [i nodes(2) dE(1,1)];
            end
        else
            G{i} = [i dE(1,1) dE(2,1)];
        end
        
        continue;
    end
    
    if (isempty(dE))
        
        D = [(1:size(Cent,1))' sqrt((Cent(:,1)-linelist(i).c(1)).^2 + (Cent(:,2)-linelist(i).c(2)).^2)];
        D = sortrows(D,2);
        
        nodes = D(:,1);
        nodes = nodes(nodes ~= i);
        
        if (isempty(dS))
            G{i} = [i nodes(1) nodes(2)];
        elseif (size(dS,1) == 1)
            if(nodes(1) ~= dS(1,1))
                G{i} = [i nodes(1) dS(1,1)];
            else
                G{i} = [i nodes(2) dS(1,1)];
            end
        else
            G{i} = [i dS(1,1) dS(2,1)];
        end
        
        continue;
    end
    
   G{i} = [i dS(1,1) dE(1,1)];
   
end
