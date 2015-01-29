function [V] = SegmentPotential(line, label, mrf)
   if(label==1)
        if(mrf.params.allowjunc == 0)
            V = 0;
        else
            V = (0.5413 - line.prob);
        end
   else
       V = (-0.4587 + line.prob);
   end
end