function [T] = FindTerminationPoints(linelist, labels, Cliques)

    k = 0;
    T = [];
    
    for i=1:size(linelist,2)
        if(labels(i) == 1)
            clqS = Cliques{linelist(i).clq(1)};
            clqS = clqS(clqS ~= i);
            
            if(~any(labels(clqS)))
                k = k+1;
                T(k,:) = [i 1];
            end
            
            clqE = Cliques{linelist(i).clq(2)};
            clqE = clqE(clqE ~= i);
            
            if(~any(labels(clqE)))
                k = k+1;
                T(k,:) = [i 2];
            end
            
        end
    end

end