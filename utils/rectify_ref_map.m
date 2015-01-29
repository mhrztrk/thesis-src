%%
for i=5:size(ref,2)
    for j=1:size(ref{i},2)
        
        if(~isempty(endp{i}))
            s_dist_ep = ((endp{i}(:,1) - ref{i}{j}.s(1)).^2 + (endp{i}(:,2) - ref{i}{j}.s(2)).^2);  
        else
            s_dist_ep = Inf;
        end
        
        if(~isempty(junc{i}))
            s_dist_jn = ((junc{i}(:,1) - ref{i}{j}.s(1)).^2 + (junc{i}(:,2) - ref{i}{j}.s(2)).^2);
        else
            s_dist_jn = Inf;
        end
        
        [min_ep, ind_ep] = min(s_dist_ep);
        [min_jn, ind_jn] = min(s_dist_jn);
        
        if(min_ep < min_jn)
            lineSet{i}{j}.s = endp{i}(ind_ep,:);
            lineSet{i}{j}.conn_status(1) = 0;
            lineSet{i}{j}.conn_through(1) = ind_ep;
        else
            lineSet{i}{j}.s = junc{i}(ind_jn,:);
            lineSet{i}{j}.conn_status(1) = 1;
            lineSet{i}{j}.conn_through(1) = ind_jn;            
        end
          
        if(~isempty(endp{i}))
            e_dist_ep = ((endp{i}(:,1) - ref{i}{j}.e(1)).^2 + (endp{i}(:,2) - ref{i}{j}.e(2)).^2);   
        else
            e_dist_ep = Inf;
        end
        
        if(~isempty(junc{i}))
            e_dist_jn = ((junc{i}(:,1) - ref{i}{j}.e(1)).^2 + (junc{i}(:,2) - ref{i}{j}.e(2)).^2); 
        else
            e_dist_jn = Inf;
        end
        
        [min_ep, ind_ep] = min(e_dist_ep);
        [min_jn, ind_jn] = min(e_dist_jn);
        
        if(min_ep < min_jn)
            lineSet{i}{j}.e = endp{i}(ind_ep,:);
            lineSet{i}{j}.conn_status(2) = 0;
            lineSet{i}{j}.conn_through(2) = ind_ep;
        else
            lineSet{i}{j}.e = junc{i}(ind_jn,:);
            lineSet{i}{j}.conn_status(2) = 1;
            lineSet{i}{j}.conn_through(2) = ind_jn;            
        end
        
    end
end

clear min_ep ind_ep min_jn ind_jn s_dist_ep s_dist_jn e_dist_ep e_dist_jn i j

%%

for i=5:8
    
    figure;imshow(img{i})
    hold on
    
    for j=1:size(ref{i},2)
        
        line([lineSet{i}{j}.s(1) lineSet{i}{j}.e(1)], [lineSet{i}{j}.s(2) lineSet{i}{j}.e(2)],'Color',[0 0 1],'LineWidth',3);
        
        if(lineSet{i}{j}.conn_status(1) == 0)
            scatter(lineSet{i}{j}.s(1), lineSet{i}{j}.s(2), 'r','.','LineWidth', 3);
        else
            scatter(lineSet{i}{j}.s(1), lineSet{i}{j}.s(2), 'g','.','LineWidth', 3);
        end
        
        if(lineSet{i}{j}.conn_status(2) == 0)
            scatter(lineSet{i}{j}.e(1), lineSet{i}{j}.e(2), 'r','.','LineWidth', 3);
        else
            scatter(lineSet{i}{j}.e(1), lineSet{i}{j}.e(2), 'g','.','LineWidth', 3);
        end
        
    end
end

clear i j
%%

for i=1:4
    DrawLabelledLineList(img{i}, linelist{i}, labels{i});  
end

%%


