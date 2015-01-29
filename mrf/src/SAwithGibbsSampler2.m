%%
function [labels_GlbMin] = SAwithGibbsSampler2(linelist, iLabels, nstep, nlevel, iTemp, mrf, Cliques)

    E_new = Inf;
    E_prv = Inf;
    E_GlbMin = Inf;
    deltaE = Inf;

    T = iTemp;
    c = 0.5;
    K = 1;
       
    labels        = iLabels;
    labels_GlbMin = iLabels; 

    nclq = size(Cliques,1);
   
    datetime = datestr(now, '(yyyy.mm.dd)_(HH.MM.SS)');
    
    while(deltaE >= 0.000 && K <= nlevel)
        
        LogFile = fopen(sprintf('MRF/mrf_optimization_%d_%s.log', K, datetime),'w+');
        
        for k=1:nstep
            
            for i=1:nclq

                ngrp = size(Cliques{i},2);
                
                labels_allp = zeros(ngrp, 2^ngrp);
                for m=1:(2^ngrp)
                    labels_allp(:,m)=de2bi((m-1),ngrp)';
                end
                
                labels_l = repmat(labels,1,(2^ngrp));
                labels_l(Cliques{i},:) = labels_allp;
                
                for m=1:(2^ngrp)
                    prob(m) = 0;
                    for n=1:ngrp
                        prob(m) = prob(m) + exp(-MRFLocalEnergy(Cliques{i}(n) , linelist, labels_l(:,m), Cliques, mrf));
                    end
                end

                prob = prob ./ sum(prob);
                sumE = 0;
                r = rand(1);
                for m=1:(2^ngrp)
                    sumE = sumE + prob(m);
                    if(r < sumE)
                        labels_gibbs = labels_l(:,m);
                        break;
                    end
                end

                % simulated annealing
                E = CalculateEnergyLL(linelist, labels_gibbs, Cliques, mrf);
                if(E < E_prv)
                    % accept change
                    E_new = E;
                    labels = labels_gibbs;
                else
                    dE = (E - E_prv);
                    % Z = exp(-E/T) + exp(-E_old/T);
                    Z = 1;
                    p = (1/Z) * exp(-dE/T);
                    r = rand(1);
                    if(r < p)
                        % accept change
                        E_new = E;
                        labels = labels_gibbs;
                    end

                end

                if(E_new < E_GlbMin)
                    E_GlbMin = E_new;
                    labels_GlbMin = labels;
                end

                deltaE = abs(E_new - E_prv);    % stop when energy change is small    
                E_prv = E_new;
                
                fwrite(LogFile, sprintf('\t__%4d__%02d(%3d), E_new = %.2f, E_prev = %.2f, E_glb = %.2f\n', i, K, k, E_new, E_prv, E_GlbMin));
                
            end
            
            fprintf('_%02d(%3d), Energy = %.3f, deltaE = %.3f\n', K, k, E_new, deltaE);
            
            fwrite(LogFile, sprintf('__%02d(%3d), E_new = %.2f, E_prev = %.2f, E_glb = %.2f\n', K, k, E_new, E_prv, E_GlbMin));   
        end
        
        labels = labels_GlbMin;
        
        T = T * c;      % Decrease temperature
        K = K + 1;      % Advance iteration counter

        fclose(LogFile);
        
    end            
end
