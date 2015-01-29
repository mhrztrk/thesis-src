
%   S: Segments 
%   L: Link labels (mx1)([0,1])
%   C: Clique Set
%   O: Observations (mx1)([0 1])
%
function E = CalculateEnergy(S, L, C, O, mrf)

    V  = SegmentPotentials(L, O, mrf);      % V(d_i|l_i)
    Vc = CliquePotentials(S, L, C, mrf);     % V_c(l)

    E = sum(V) + sum(Vc);
    
end