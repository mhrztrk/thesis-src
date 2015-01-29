function [V] = SegmentPotentialsLL(linelist, labels, mrf)

%     V = zeros(size(labels));
%     
%     for i=1:size(labels,1)
%        V(i) = SegmentPotential(linelist(i), labels(i), mrf);
%     end

    
   V = 0;
   for i=1:size(labels,1)
       if(labels(i)==1)
            if(mrf.params.allowjunc == 0)
              
            else
                V = V + linelist(i).SegPot(2);
            end
       else
           V = V + linelist(i).SegPot(1);
       end
   end
  
    
   
end