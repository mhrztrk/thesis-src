function [LabelMinEx, LineSetEx, CliqSetEx] = ConnectTerminationPoints(LineSet, CliqSet, LabelMin, mrf, centroids, pmap)

    LineSetEx = LineSet;
    CliqSetEx = CliqSet;
    LabelMinEx = LabelMin;    

    mrf.params.allowjunc = 1;

    extConn = [];
    nconn = 0;
    mxDist = 45;
    
    for i=1:size(CliqSet,1)

        labels_tmp = LabelMin;
        plist = CliqSet{i}((LabelMin(CliqSet{i})==1));

        if(size(plist,2) == 1)
            
            nlist = CliqSet{i}((LabelMin(CliqSet{i})==0));

            if(LineSetEx(plist(1)).clq(1) == i)
                term_pt = LineSetEx(plist(1)).s;
                term_cnod = LineSetEx(plist(1)).adj(1);
            else
                term_pt = LineSetEx(plist(1)).e;
                term_cnod = LineSetEx(plist(1)).adj(2);
            end
            
            nodes = FindAllCloseNodes(plist(1), term_pt, LineSetEx, mxDist, centroids);
            
            if(isempty(nodes))
                continue
            end
            
            line = [];
            exist = [];
            for j=1:size(nodes,1)
                [exist(j) line(j)] = LineExist(LineSetEx, term_pt, centroids(nodes(j),:));
            end
         
            for j=1:size(nodes,1)
                if(exist(j)==0)
                    % Add new linesDr
                    [LineSetEx CliqSetEx linenr] = AddNewLine(LineSetEx, CliqSetEx, term_pt, centroids(nodes(j),:), term_cnod, nodes(j), pmap);
                    exist(j) = 1;
                    line(j) = linenr;
                    LabelMinEx(linenr) = 0;
                end
            end
            
            Er = [];
            Ei = CalculateEnergyLL(LineSetEx, LabelMinEx, CliqSetEx, mrf);
            for j=1:size(nodes,1)
                LabelMinEx(line(j)) = 1;
                Er(j) = CalculateEnergyLL(LineSetEx, LabelMinEx, CliqSetEx, mrf);
                LabelMinEx(line(j)) = 0;
            end
            
            [val ind] = min(Er);
            if(val < Ei)
                LabelMinEx(line(ind)) = 1;
                %DrawLabelledLineList(pmap, LineSetEx, LabelMinEx);
            end
%%            
%             if(size(nlist,2) > 1) 
%             
%                 combos = combntns(1:size(nlist,2), 2);
%                 pts = zeros(size(nlist,2),2);
% 
%                 for k=1:size(nlist,2)
%                     if(LineSet(nlist(k)).clq(1) == i)
%                         pts(k,:) = LineSet(nlist(k)).e;
%                     else
%                         pts(k,:) = LineSet(nlist(k)).s;
%                     end
%                 end
% 
%                 for k=1:size(combos,1)
%                     % check for intersection
%                     [x y] = FindIntersectionPoint(LineSet(plist(1)).s, LineSet(plist(1)).e, pts(combos(k,1),:),pts(combos(k,2),:));
%                     
%                     len = sqrt((pts(combos(k,1),1)-pts(combos(k,2),1))^2+(pts(combos(k,1),2)-pts(combos(k,2),2))^2);
%                     d11 = sqrt((x-pts(combos(k,1),1))^2+(y-pts(combos(k,1),2))^2);
%                     d12 = sqrt((x-pts(combos(k,2),1))^2+(y-pts(combos(k,2),2))^2);
%             
%                     if((d11+d12) <= len)
%                         
%                         nconn = nconn + 1;
%                         ConnS(nconn,:) = term_pt;
%                         ConnE(nconn,:) = [x y];
%                         % disp('found');
%                     end
%                     
%                     
%                     
%                 end
%             end
%             
%             Er = zeros(size(nlist));
%             for j=1:size(nlist,2)
%                 labels_tmp(nlist(j)) = 1;
%                 Er(j) = CalculateEnergyLL(LineSet, labels_tmp, CliqSet, mrf);
%                 labels_tmp(nlist(j)) = 0;
%             end
% 
%             [val ind] = min(Er);
%             if(val < Ef)
%                 extConn = [extConn nlist(ind)];
%             end
%%
        end

    end

    labelsExt = LabelMin;
    labelsExt(extConn) = ones(size(extConn)); 

end


function [x y] = FindIntersectionPoint(sp1,ep1,sp2,ep2)

    m1 = (sp1(2)-ep1(2))/(sp1(1)-ep1(1));
    m2 = (sp2(2)-ep2(2))/(sp2(1)-ep2(1));
   
    c1 = sp1(2) - m1*sp1(1);
    c2 = sp2(2) - m2*sp2(1);
    
    
    if(m1 == m2) % parallel lines
        x = [];
        y = [];
    elseif(isinf(m1))
        x = sp1(1);
        y = c2 + m2*x;
    elseif(isinf(m2)) 
        x = sp2(1);
        y = c1 + m1*x;
    else 
        x = (c1 - c2)/(m2 - m1);
        y = c1 + m1*x;
    end
    
end


function [nodes] = FindAllCloseNodes(seg, term_pt, LineSet, mxDist, centroids)

   D = [(1:size(centroids,1))' sqrt((centroids(:,1)-term_pt(1)).^2 + (centroids(:,2)-term_pt(2)).^2)];
   D = sortrows(D,2);
   nodes = D(:,1);
   
   nodes = nodes(2:sum(D(:,2) <= mxDist));
   
   if(all(LineSet(seg).s == term_pt))
       ep = LineSet(seg).e;
   else
       ep = LineSet(seg).s;
   end
   
   ang = GetAngle((ep(1) - term_pt(1)),(ep(2) - term_pt(2)));
   
   for i=1:size(nodes,1)
       connAng(i) = GetAngle((centroids(nodes(i),1) - term_pt(1)),(centroids(nodes(i),2) - term_pt(2)));
   end
   
   defl = abs(connAng - ang);
   defl(defl>180) = 360 - defl(defl>180);
   
   nodes = nodes(defl>90);
   
end


function [exist line] = LineExist(LineSet, pt1, pt2)

    for i=1:size(LineSet,2)
    
        if((all(LineSet(i).s == pt1) && all(LineSet(i).e == pt2)) || ...
               (all(LineSet(i).e == pt1) && all(LineSet(i).s == pt2)))
            exist = 1;
            line = i;
            break;
        else
            exist = 0;
            line = -1;
        end
    end
end


function [LineSetEx CliqSetEx nconn] = AddNewLine(LineSet, CliqSet, sp, ep, sCC, eCC, pmap)

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

