
nconn = size(linelist, 2);

tmp = [linelist.c];

Cent(:,1) = tmp(1:2:end);
Cent(:,2) = tmp(2:2:end);

for i=1:nconn

   D = [(1:size(Cent,1))' sqrt((Cent(:,1)-linelist(i).c(1)).^2 + (Cent(:,2)-linelist(i).c(2)).^2)];
    
   D = sortrows(D,2);
    
   nodes = D(:,1);
   
   G{i} = nodes(1:5);
   
end

G = G';