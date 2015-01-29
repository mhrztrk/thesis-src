function [labels] = MRFRelaxWithICM(linelist, iLabels, nstep, mrf, Cliques)
      
    mrf.CliqPotMin = -1 * (-mrf.params.Kl*2 + mrf.params.Kc*sind(180));
    
    for i=1:size(linelist,2)
        
        linelist(i).SegPot(1) = (-0.4587 + linelist(i).prob);
        
        if(mrf.params.allowjunc == 0)
            linelist(i).SegPot(2) = 0;
        else
            linelist(i).SegPot(2) = (0.5413 - linelist(i).prob);
        end
        
    end

    labels        = iLabels;

    ngrp = 8;
    [G] = CreateGroupMatrix(linelist, ngrp);

    labels_allp = zeros(ngrp, 2^ngrp);
    for m=1:(2^ngrp)
        labels_allp(:,m)=de2bi((m-1),ngrp)';
    end
    
    E = CalculateEnergyLL(linelist, labels, Cliques, mrf);
    fprintf('Initial Energy = %.3f\n', E);
    
    for k=1:nstep

        for i=1:size(G,1)

            labels_l = repmat(labels,1,(2^ngrp));
            labels_l(G{i},:) = labels_allp;

            for m=1:(2^ngrp)
                E_loc(m) = 0;
                for n=1:ngrp
                    E_loc(m) = E_loc(m) + MRFLocalEnergy(G{i}(n) , linelist, labels_l(:,m), Cliques, mrf);
                end
            end

            [~, ind] = min(E_loc);
            labels = labels_l(:,ind);

        end

        E = CalculateEnergyLL(linelist, labels, Cliques, mrf);
        fprintf('_%02d, Energy = %.3f\n', k, E);

    end

end


function [Grp] = CreateGroupMatrix(linelist, nset)

    nconn = size(linelist, 2);

    tmp = [linelist.c];

    Cent(:,1) = tmp(1:2:end);
    Cent(:,2) = tmp(2:2:end);

    Grp = cell(nconn,1);
    
    for i=1:nconn
       D = [(1:size(Cent,1))' sqrt((Cent(:,1)-linelist(i).c(1)).^2 + (Cent(:,2)-linelist(i).c(2)).^2)];
       D = sortrows(D,2);
       nodes = D(:,1);
       Grp{i} = nodes(1:nset);
    end
    
end
