function [V] = SegmentPotentials(labels, D, mrf)

    V = zeros(size(labels));
    
    for i=1:size(labels,1)
       if(labels(i)==0)
           V(i) = 0;
       else
           if(D(i) < mrf.params.t1)
               V(i) = 0;
           elseif (D(i) >= mrf.params.t1 && D(i) <= mrf.params.t2)
               V(i) = (D(i) - mrf.params.t1)/(mrf.params.t2 - mrf.params.t1);
           else
               V(i) = 1;
           end
       end
    end
    
end