function [linelist, Cliques, CC] = KovesiLineModel2MRFNetwork(rdclseglist, pmap, mxDist, minProb)

    k = 0;
    CC = [];
    CC_isendp = [];
    
    for j=1:size(rdclseglist,2)
        for m=1:size(rdclseglist{j},1)
            
            k = k + 1;
            
            CC(k,:) = rdclseglist{j}(m,:);
            
            if(m == 1)
                ang = GetAngle((CC(m,1)-rdclseglist{j}(m+1,1)), (CC(k,2)-rdclseglist{j}(m+1,2)));
                CC_isendp(k,:) = [1 ang k+1 0];
            elseif(m == size(rdclseglist{j},1))
                ang = GetAngle((CC(m,1)-rdclseglist{j}(m-1,1)), (CC(k,2)-rdclseglist{j}(m-1,2)));
                CC_isendp(k,:) = [1 ang k-1 0];
            else
                CC_isendp(k,:) = [0 0 k-1 k+1];
            end
        end
    end
    
    [linelist, Cliques] = SORMNetworkConstructEpOnly(CC, CC_isendp, pmap, mxDist, minProb);
     
end


function [linelist, Cliques] = SORMNetworkConstructEpOnly(CC, CC_isendp, pmap, mxDist, minProb)

    ang_th  = 16;
    
    Adj = zeros(size(CC,1));
    
    isendp = CC_isendp(:,1);
    endpang = CC_isendp(:,2);
    connpt = CC_isendp(:,3:4);
    
    for i=1:size(CC,1)
                
        % get points close to point i
        D = sqrt((CC(:,1)-CC(i,1)).^2 + (CC(:,2)-CC(i,2)).^2);

        clsPts = find(D < mxDist);

        % remove self
        clsPts = clsPts(clsPts ~= i); 
        
        if(~any(connpt(i,1)== clsPts))
            clsPts = [clsPts;connpt(i,1)];
        end
        
        if(isendp(i) == 0)
            if(~any(connpt(i,2)== clsPts))
                clsPts = [clsPts;connpt(i,2)];
            end
        end
        
        prob = zeros(1,size(clsPts,1));
        Ang  = zeros(1,size(clsPts,1)); 

        for j=1:size(clsPts,1)
           Ang(j) = GetAngle((CC(i,1)-CC(clsPts(j),1)), (CC(i,2)-CC(clsPts(j),2))); % link orientation    
           prob(j) = GetLineProb(pmap, CC(i,:), CC(clsPts(j),:));
        end

        orien_diff = abs(Ang - endpang(i));
        orien_diff(orien_diff>180) = 360 - orien_diff(orien_diff>180); 

        Ang = Ang(((prob > minProb) & (orien_diff > 90)) | ((clsPts ~= connpt(i,1))' | (clsPts ~= connpt(i,2))'));
        clsPts = clsPts(((prob > minProb) & (orien_diff > 90)) | ((clsPts ~= connpt(i,1))' | (clsPts ~= connpt(i,2))'));
        
        % if there are more than one point on the same orientation, take the
        % closest one only

        if(~isempty(clsPts))

            Mat = sortrows([clsPts D(clsPts) Ang'],2);

            k = 1;

            AllocReg = [];
            AllocReg(1,:) = [(Mat(1,3)-(ang_th/2)) (Mat(1,3)+(ang_th/2))];
            clsPtsTmp = Mat(1,1);

            for m=2:size(clsPts,1)
                notA = 0;
                for n=1:k
                    if(Mat(m,3) > AllocReg(n,1) && Mat(m,3) < AllocReg(n,2))
                        notA = 1;
                    end
                end
                if(notA == 0)
                    k = k + 1;
                    AllocReg(k,:) = [(Mat(m,3)-(ang_th/2)) (Mat(m,3)+(ang_th/2))];
                    clsPtsTmp(k) = Mat(m,1);
                end
            end

            clsPts = clsPtsTmp';

            for j=1:size(clsPts,1)
               Adj(i,clsPts(j)) = 1;
               Adj(clsPts(j),i) = 1;          
            end
        end
    end
    
    
    [X Y] = find(triu(Adj)==1); 
  
    for i=1:size(X,1)
        
        linelist(i).s = CC(X(i),:); % link start point
        linelist(i).e = CC(Y(i),:); % link end point
        
        linelist(i).c = (linelist(i).s + linelist(i).e)/2;
        
        linelist(i).len = sqrt((linelist(i).s(1) - linelist(i).e(1))^2 + ...
                               (linelist(i).s(2) - linelist(i).e(2))^2);    % link length

                           
        linelist(i).ang = GetAngle((linelist(i).e(1) - linelist(i).s(1)), ...
                            (linelist(i).e(2) - linelist(i).s(2))); % link orientation
                        
        linelist(i).adj = [X(i) Y(i)];
    
        linelist(i).nclq = 2; 
        
        linelist(i).prob = GetLineProb(pmap, linelist(i).s, linelist(i).e);
        
        linelist(i).nr = i;
    end
   
    maxLen = max([linelist.len]);
    for i=1:size(X,1)
        linelist(i).len = linelist(i).len / maxLen;
    end
    
    [Cliques,linelist] = FindConnCliques(linelist, Adj);
   
    fileID = fopen('SormNodes.txt', 'w+');
    
    for j=1:size(linelist,2)
        fprintf(fileID, '%4d - ( %5.1f , %5.1f )\n', j, linelist(j).c(1), linelist(j).c(2));
    end
    
end

%
% figure;imshow(clres);hold on;
% 
% scatter(cc(:,2),cc(:,1),'r','x');
% 
% for i=1:size(adj,1)
%     for j=1:size(adj,1)
%         if(adj(i,j) == 1)
%             line([cc(i,2) cc(j,2)],[cc(i,1) cc(j,1)],'linewidth',2);    
%         end
%     end
% end
%


function [Clqs, linelist] = FindConnCliques(linelist, Adj)

    Adj = tril(Adj);

    for i=1:size(linelist,2)
        linelist(i).ClqS = 0;
        linelist(i).ClqE = 0;
    end
    
    tmp = [linelist.adj];
    AdjV(:,1) = tmp(1:2:end);
    AdjV(:,2) = tmp(2:2:end);
    
    k = 0;
    
    for i=1:size(linelist,2)
        
        connS  = [find((AdjV(:,1) == AdjV(i,1)));find((AdjV(:,2) == AdjV(i,1)))]; connS = connS(connS~=i);
        connE  = [find((AdjV(:,1) == AdjV(i,2)));find((AdjV(:,2) == AdjV(i,2)))]; connE = connE(connE~=i);
        
        if(linelist(i).ClqS == 0)
            
            k = k + 1;
            Clqs{k} = i;
            if(~isempty(connS))
                for j=1:size(connS,1)
                    if(all(linelist(i).s == linelist(connS(j)).s))
                        if(linelist(connS(j)).ClqS == 0)
                            linelist(connS(j)).ClqS = 1;
                            linelist(connS(j)).clq(1) = k;
                            Clqs{k} = [Clqs{k} connS(j)];
                        end
                    elseif(all(linelist(i).s == linelist(connS(j)).e))
                        if(linelist(connS(j)).ClqE == 0)
                            linelist(connS(j)).ClqE = 1;
                            linelist(connS(j)).clq(2) = k;
                            Clqs{k} = [Clqs{k} connS(j)];
                        end
                    else
                        error('wtf!!!');
                    end
                        
                end
            end
            linelist(i).ClqS = 1;
            linelist(i).clq(1) = k;
        end
        
        if(linelist(i).ClqE == 0)
            
            k = k + 1;
            Clqs{k} = i;
            if(~isempty(connE))
                for j=1:size(connE,1)
                    if(all(linelist(i).e == linelist(connE(j)).s))
                        if(linelist(connE(j)).ClqS == 0)
                            linelist(connE(j)).ClqS = 1;
                            linelist(connE(j)).clq(1) = k;
                            Clqs{k} = [Clqs{k} connE(j)];
                        end
                    elseif(all(linelist(i).e == linelist(connE(j)).e))
                        if(linelist(connE(j)).ClqE == 0)
                            linelist(connE(j)).ClqE = 1;
                            linelist(connE(j)).clq(2) = k;
                            Clqs{k} = [Clqs{k} connE(j)];
                        end
                    else
                        error('wtf!!!');
                    end
                        
                end
            end
            linelist(i).ClqE = 1;
            linelist(i).clq(2) = k;
        end 
    end

    Clqs = Clqs';
    
end


