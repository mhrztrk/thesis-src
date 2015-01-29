function E = CalculateNetworkEnergy(linelist, labels, clique, mrf)
    
    mrf.CliqPotMin = -1 * (-mrf.params.Kl*2 + mrf.params.Kc*sind(180));
    
    for i=1:size(linelist,2)
        
        linelist(i).SegPot(1) = (-0.4587 + linelist(i).prob);
        
        if(mrf.params.allowjunc == 0)
            linelist(i).SegPot(2) = 0;
        else
            linelist(i).SegPot(2) = (0.5413 - linelist(i).prob);
        end
        
    end
   
     
    Vtot  = mrf.params.eratio * SegmentPotentialsLL(linelist, labels, mrf);      % V(d_i|l_i)
    
    Vc = CliquePotentialsLL(linelist, labels, clique, mrf);     % V_c(l)
  
    E = Vtot + sum(Vc);
    
end