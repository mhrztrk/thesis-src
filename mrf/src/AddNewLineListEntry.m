function [LineSetEx CliqSetEx nconn] = AddNewLineListEntry(LineSet, CliqSet, sp, ep, sCC, eCC, pmap)

    LineSetEx = LineSet;
    CliqSetEx = CliqSet;
    
    nconn = size(LineSet, 2);
    nconn = nconn + 1;
    
    LineSetEx(nconn).s = sp; % link start point
    LineSetEx(nconn).e = ep; % link end point

    LineSetEx(nconn).c = (LineSetEx(nconn).s + LineSetEx(nconn).e)/2;

    LineSetEx(nconn).len = sqrt((LineSetEx(nconn).s(1) - LineSetEx(nconn).e(1))^2 + ...
                           (LineSetEx(nconn).s(2) - LineSetEx(nconn).e(2))^2);    % link length


    LineSetEx(nconn).ang = GetAngle((LineSetEx(nconn).e(1) - LineSetEx(nconn).s(1)), ...
                        (LineSetEx(nconn).e(2) - LineSetEx(nconn).s(2))); % link orientation

    LineSetEx(nconn).adj = [sCC eCC];

    LineSetEx(nconn).nclq = 2; 

    LineSetEx(nconn).prob = GetLineProb(pmap, sp, ep);

    LineSetEx(nconn).nr = nconn;

    for i=1:size(CliqSet,1)
        for j=1:size(CliqSet{i},2)
            
            if(LineSet(CliqSet{i}(j)).clq(1)==i)
                pt = LineSet(CliqSet{i}(j)).s;
            else
                pt = LineSet(CliqSet{i}(j)).e; 
            end
            
            if(all(pt == sp))
                LineSetEx(nconn).ClqS = 1;
                LineSetEx(nconn).clq(1) = i;
                CliqSetEx{i}(size(CliqSet{i},2)+1) = nconn;
                break;
            end
        end 
    end
    
    for i=1:size(CliqSet,1)
        for j=1:size(CliqSet{i},2)
            
            if(LineSet(CliqSet{i}(j)).clq(1)==i)
                pt = LineSet(CliqSet{i}(j)).s;
            else
                pt = LineSet(CliqSet{i}(j)).e; 
            end
            
            if(all(pt == ep))
                LineSetEx(nconn).ClqSE = 1;
                LineSetEx(nconn).clq(2) = i;
                CliqSetEx{i}(size(CliqSet{i},2)+1) = nconn;
                break;
            end
        end 
    end
    
    [LineSetEx] = UpdateLineLength(LineSetEx);
      
end

function [LineSet] = UpdateLineLength(LineSet)

    for i=1:size(LineSet,2)
        
        LineSet(i).len = sqrt((LineSet(i).s(1) - LineSet(i).e(1))^2 + ...
                               (LineSet(i).s(2) - LineSet(i).e(2))^2);    % link length


    end
   
    maxLen = max([LineSet.len]);
    for i=1:size(LineSet,2)
        LineSet(i).len = LineSet(i).len / maxLen;
    end
    
end
