labelsOpt = labels_GlbMin_2;
labelsMin = labelsOpt;

[T] = FindTerminationPoints(linelist, labelsOpt, Cliques);

mrf.params.allowjunc = 1;

nline = size(linelist,2);
nclq  = size(Cliques, 1);

tmp = [linelist.c];
C = [tmp(1:2:end)' tmp(2:2:end)'];

linelist2 = linelist;
Cliques2 = Cliques;

figure;imshow(img(:,:,1:3),[]); hold on;

for i=1:size(T,1)
    
    ind = T(i,1);

    if(T(i,2)==1)
        pt = linelist(ind).s;
        clq = linelist(ind).clq(1);
    else
        pt = linelist(ind).e;
        clq = linelist(ind).clq(2);
    end
    
    conns = Cliques{clq};
    conns = conns(conns~=ind);
    
    for j=1:(size(conns,2)-1)
        
        if(linelist(conns(j)).clq(1)~=clq)
            clq_j = linelist(conns(j)).clq(1);
        else
            clq_j = linelist(conns(j)).clq(2);
        end
        
        for k=(j+1):size(conns,2)
   
            if(linelist(conns(k)).clq(1)~=clq)
                clq_k = linelist(conns(k)).clq(1);
            else
                clq_k = linelist(conns(k)).clq(2);
            end
            
            comm = intersect(Cliques{clq_k},Cliques{clq_j});
            
            if(labelsOpt(comm) == 1)
                
                line([linelist(ind).s(2) linelist(ind).e(2)], [linelist(ind).s(1) linelist(ind).e(1)],'Color',[0 0 1],'LineWidth',5);
                
                line([linelist(conns(j)).s(2) linelist(conns(j)).e(2)], [linelist(conns(j)).s(1) linelist(conns(j)).e(1)],'Color',[0 1 0],'LineWidth',5);
                line([linelist(conns(k)).s(2) linelist(conns(k)).e(2)], [linelist(conns(k)).s(1) linelist(conns(k)).e(1)],'Color',[0 1 0],'LineWidth',5);
                
                line([linelist(comm).s(2) linelist(comm).e(2)], [linelist(comm).s(1) linelist(comm).e(1)],'Color',[1 0 0],'LineWidth',5);
                
                m1 = (linelist(ind).s(1) - linelist(ind).e(1))/(linelist(ind).s(2) - linelist(ind).e(2));
                m2 = (linelist(comm).s(1) - linelist(comm).e(1))/(linelist(comm).s(2) - linelist(comm).e(2));
                
                s1 = linelist(ind).s;
                s2 = linelist(comm).s;
                
                y = (m1*s1(2) - m2*s2(2) - (s1(1) - s2(1)))/(m1 - m2);
                x = m1*(y - s1(2)) + s1(1);
                
                dist2s = sqrt((linelist(comm).s(1)-x)^2 + (linelist(comm).s(2)-y)^2);
                dist2e = sqrt((linelist(comm).e(1)-x)^2 + (linelist(comm).e(2)-y)^2);
                len = sqrt((linelist(comm).e(1)-linelist(comm).s(1))^2 + (linelist(comm).e(2)-linelist(comm).s(2))^2);
                
                
                if(max(dist2e,dist2s) < len)
              
                    linelist2(nline + 1).nr = (nline + 1);
                    linelist2(nline + 1).s = pt;
                    linelist2(nline + 1).e = [x y];
                    linelist2(nline + 1).c = (linelist2(nline + 1).s + linelist2(nline + 1).e)/2;
                    
                    linelist2(nline + 2).nr = (nline + 2);
                    linelist2(nline + 2).s = linelist(comm).s;
                    linelist2(nline + 2).e = [x y];
                    linelist2(nline + 2).c = (linelist2(nline + 2).s + linelist2(nline + 2).e)/2;
                    
                    linelist2(comm).s = [x y]; 
                    linelist2(comm).e = linelist(comm).e;
                    linelist2(comm).c = (linelist(comm).s + linelist(comm).e)/2;
                    
                    Cliques2{nclq+1} = [(nline+1) (nline+2) comm];
                    
                    Cliques2{linelist(comm).clq(1)} = Cliques2{linelist(comm).clq(1)}(Cliques2{linelist(comm).clq(1)} ~= comm); 
                    
                    Cliques2{clq} = [Cliques2{clq} (nline+1)];
                    
                    scatter(y, x, 'o', 'b', 'LineWidth',5);
                    [y x]
                end
                
            end
            
        end
    end
   
    %%
    Ei = CalculateEnergyLL(linelist, labelsOpt, Cliques, mrf);
    E = [];
    labelsTmp = labelsOpt;

    for j=1:size(conns,2)
        labelsTmp(conns(j)) = 1;
        E(j) = CalculateEnergyLL(linelist, labelsTmp, Cliques, mrf);
        labelsTmp(conns(j)) = 0;
    end
    
    [Emin, c] = min(E);
    if(Emin < Ei)
        labelsMin(conns(c)) = 1;
    end
    %%
    
    %%
    clq_p = linelist(ind).clq(linelist(ind).clq ~= clq); 
    
    ind_p = Cliques{clq_p}((labelsOpt(Cliques{clq_p}) == 1) & (Cliques{clq_p} ~= ind)');
    
    if(isempty(ind_p) || (max(size(ind_p)) > 1))
        continue;
    end
    
    E = [];
    labelsTmp = labelsOpt;
    
    conns_p = Cliques{clq_p}((Cliques{clq_p} ~= ind_p) & (Cliques{clq_p} ~= ind)); 
    
    labelsTmp(ind) = 0;
    
    for j=1:size(conns_p,2)
        labelsTmp(conns_p(j)) = 1;
        E(j) = CalculateEnergyLL(linelist, labelsTmp, Cliques, mrf);
        labelsTmp(conns_p(j)) = 0;
    end
    
    [Emin_p, c] = min(E);
    if(Emin_p < min(Emin, Ei))
        labelsMin(ind) = 0;
        labelsMin(conns_p(c)) = 1;
    end
    
    %%
    
    
%     nodes = 1:nline;
%     nodes = nodes(labelsOpt==1);
%     nodes = nodes(nodes~=ind);
%     
%     D = [nodes' sqrt((C(nodes,1) - pt(1)).^2 + (C(nodes,2) - pt(2)).^2)]; 
%     D = sortrows(D,2);
%     D = D(D(:,2)<40,:);

end

