function E = MRFLocalEnergy(node, linelist, labels, Cliques, mrf)

        %V  = SegmentPotential(linelist(node), labels(node), mrf);
        
        if(labels(node) == 0)
            V  = linelist(node).SegPot(1);
        else
            V  = linelist(node).SegPot(1);
        end

        if(mrf.sorm == 1)
            Vc =  0.5*(SormCliquePotential(linelist, labels, Cliques{linelist(node).clq(1)}, mrf) + ...
                        SormCliquePotential(linelist, labels, Cliques{linelist(node).clq(2)}, mrf));            
        else            
            if(linelist(node).nclq == 2) 
                Vc =  0.5*(CliquePotential(linelist, labels, Cliques{linelist(node).clq(1)}, mrf) + ...
                        CliquePotential(linelist, labels, Cliques{linelist(node).clq(2)}, mrf));
            else
                Vc =  CliquePotential(linelist, labels, Cliques{linelist(node).clq(1)}, mrf);
            end
        end

        E = (Vc + mrf.params.eratio* V);
end