function [V] = CliquePotentials(Segments, labels, Cliques, mrf)
    
    V = zeros(size(Cliques));
    
    for c=1:size(Cliques,1)
        
        ind = find(labels(Cliques{c})==1);
        
        % ?i?c , l_i=0
        if(size(ind,1)==0)
            V(c) = 0;
            
        % ?!i?c | l_i=1
        elseif(size(ind,1)==1) % only one segment labelled as road.
            p = Cliques{c}(ind(1));
            V(c) = mrf.params.Ke - mrf.params.Kl*Segments.len(p);
            
        % ?!(i,j)?c^2 | (l_i=l_j=1 , R_ij>?/2)
        elseif(size(ind,1)==2) % only two segment labelled as road.
            p1 = Cliques{c}(ind(1));
            p2 = Cliques{c}(ind(2));
            defl = abs(Segments.ang(p1)-Segments.ang(p2));
            if(defl <= 90)
                V(c) = -mrf.params.Kl*(Segments.len(p1)+Segments.len(p2)) + ...
                            mrf.params.Kc*sind(defl);
            else
                V(c) = mrf.params.Ki*(Segments.len(p1)+Segments.len(p2)); 
            end
            
        % in all other cases
        else
            V(c) = mrf.params.Ki*sum(Segments.len(Cliques{c}(ind)));
        end
    end
end