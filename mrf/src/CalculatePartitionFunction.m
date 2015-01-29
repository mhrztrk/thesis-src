function [Z] = CalculatePartitionFunction(Cliques, mrf, linelist,T)

    

    T = 0;
    for c=1:size(Cliques)
        T = T + 2^size(Cliques{c},1);
    end
    
    labels = zeros(1, size(linelist,2));
    
    Z = 1;
    for c=1:size(Cliques,1)
        if(mrf.sorm == 1)
            
            nseg = size(Cliques{c},2);
            
            Zt = 0;
            for i=1:2^nseg
                labels(Cliques{c}) = de2bi(i-1, nseg);
                Zt = Zt + exp( -1 * SormCliquePotential(linelist, labels, Cliques{c}, mrf) / T);
            end
            
            Z = Z * Zt*(T-2^nseg);
            
        else
            error('sorm only!!!')
        end
    end
    
end