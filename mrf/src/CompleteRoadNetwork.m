function [LineList] = CompleteRoadNetwork(LineList)

    % find lines which have unconnected endpoint
    endLines = LineList([LineList.nconn] == 1);

    for i=1:size(endLines,2)
        
        for j=1:size(LineList,2)
            
            defl1(j) = 0;
            defl2(j) = 0;
            dist(j)  = Inf;
            
            if(endLines(i).nr == LineList(j).nr)
                continue;
            end
            
            [line_i line_j nline state] = LineMinDist(endLines(i), LineList(j));
            if(state == 0)      % no connection
            elseif(state == 1)  % Li & Lj directly connected
            elseif(state == 2)  % Li divided into two line
            elseif(state == 3)  % Lj divided into two line
                
                LineList2 = ApplyLineCuttings(LineList);
                
            elseif(state == 4)  % Li & Lj connected through another line.
                if(line_i.connpt(end) ~= line_i.connpt(1))
                    defl1(j) = defl_angle(line_i, nline);
                    defl2(j) = defl_angle(nline, line_j);
                    dist(j)  = sqrt((nline.s(1)-nline.e(1))^2 + (nline.s(2)-nline.e(2))^2);
                end
            end             
        end
        
        if(~isinf(min(dist)))
            [~,ind] = min(dist);
            [line_i,line_j,~,~] = LineMinDist(endLines(i), LineList(ind));
            LineList = AddNewLine(LineList, line_i.nr, line_j.nr, line_i.connpt(end), line_j.connpt(end));
        end
        
    end
    
end
