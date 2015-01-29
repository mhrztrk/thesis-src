
%   linelist: list of all detected and derived lines  
%   labels: Link labels (mx1)([0,1])
%   clique: Clique Set
%
function E = CalculateEnergyLL(linelist, labels, clique, mrf)

    Vtot  = mrf.params.eratio * SegmentPotentialsLL(linelist, labels, mrf);      % V(d_i|l_i)
    Vc = CliquePotentialsLL(linelist, labels, clique, mrf);     % V_c(l)
  
    E = Vtot + sum(Vc);
    
end