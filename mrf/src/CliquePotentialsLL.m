function [V] = CliquePotentialsLL(linelist, labels, Cliques, mrf)
    
    V = zeros(size(Cliques));
    for c=1:size(Cliques,1)
        if(mrf.sorm == 1)
            V(c) = SormCliquePotential(linelist, labels, Cliques{c}, mrf);
        else
            V(c) = CliquePotential(linelist, labels, Cliques{c}, mrf);
        end
    end
    
end



