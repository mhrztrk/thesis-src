for i=1:15
    figure('name',sprintf('%s', img_name{i}));
    subplot(1,2,1), imshow(clres{i},[]);
    subplot(1,2,2), imshow(clres_filt{i})
end
%%

for i=1:15
    
    img_skel = bwmorph(clres_filt_th{i},'skel',Inf);
    img_bp = bwmorph(img_skel,'branchpoints');
    img_bp = bwmorph(img_bp,'dilate');
    img_skel(img_bp==1) = 0;
    img_skel = cc_threshold(img_skel, 15, 0);

    clres_filt_th_ref{i} = bwmorph(img_skel,'dilate'); 

    figure;imshow(clres_filt_th_ref{i},[])
    
    % Kovesi - EdgeLinking Toolbox
    % Link edge pixels together into lists of sequential edge points, one
    % list for each edge contour. A contour/edgelist starts/stops at an 
    % ending or a junction with another contour/edgelist.
    % Here we discard contours less than 10 pixels long.
    [edgelist, labelededgeim] = edgelink(img_skel, 10);

    % Fit line segments to the edgelists
    tol = 2;  % Line segments are fitted with maximum deviation from
              % original edge of 2 pixels.

    % Road centerline segment list          
    rdclseglist{i} = lineseg(edgelist, tol); 

    % Display the edgelists with random colours for each distinct edge 
    % in figure 2
    h = figure;
    imshow(img{i}(:,:,1:3),[]);
    drawedgelist(rdclseglist{i}, [size(img{i},1) size(img{i},2)], 1, 'rand', h); axis off 
    
end

clear img_bp i edgelist tol img_skel labelededgeim h
%%

for i=1:1
   [CC{i}, CM{i}] = sorm(clres_filt_th_ref{i}, 25, 12, 5); 
end


    
%%    
    % deleting an element from a cell array
%     rdclseglist{12} = [];
%     rdclseglist(cellfun(@(rdclseglist) isempty(rdclseglist),rdclseglist))=[];
    
for i=1:15
    
    k = 0;
    n_ep = 0;
    CC_kov{i} = [];
    CC_isendp{i} = [];
    
    for j=1:size(rdclseglist{i},2)
        for m=1:size(rdclseglist{i}{j},1)
            
            k = k + 1;
            
            CC_kov{i}(k,:) = rdclseglist{i}{j}(m,:);
            
            if(m == 1)
                ang = GetAngle((CC_kov{i}(m,1)-rdclseglist{i}{j}(m+1,1)), (CC_kov{i}(k,2)-rdclseglist{i}{j}(m+1,2)));
                CC_isendp{i}(k,:) = [1 ang];
            elseif(m == size(rdclseglist{i}{j},1))
                ang = GetAngle((CC_kov{i}(m,1)-rdclseglist{i}{j}(m-1,1)), (CC_kov{i}(k,2)-rdclseglist{i}{j}(m-1,2)));
                CC_isendp{i}(k,:) = [1 ang];
            else
                CC_isendp{i}(k,:) = [0 0];
            end
        end
    end
    
    [linelist{i}, Cliques{i}] = SORMNetworkConstructEpOnly(CC_kov{i}, CC_isendp{i}, clres{i}, 75, 0);

end

%%

%%  
for i=1:15
    
    nline = size(linelist{i},2);
    linelist_up{i} = linelist{i};
    
    for j=1:size(rdclseglist{i},2)
        for m=1:(size(rdclseglist{i}{j},1)-1)
            line_found = 0;
            for n=1:size(linelist{i},2)
                if((all(linelist{i}(n).s == rdclseglist{i}{j}(m,:)) && all(linelist{i}(n).e == rdclseglist{i}{j}(m+1,:))) || ...
                   (all(linelist{i}(n).e == rdclseglist{i}{j}(m,:)) && all(linelist{i}(n).s == rdclseglist{i}{j}(m+1,:))))
                    line_found = 1;
                end
            end
            if(line_found == 0)
                
                display('new line!')
                
                nline = nline + 1;
                
                linelist_up{i}(nline).s = rdclseglist{i}{j}(m,:); % link start point
                linelist_up{i}(nline).e = rdclseglist{i}{j}(m+1,:); % link end point

                linelist_up{i}(nline).c = (linelist_up{i}(nline).s + linelist_up{i}(nline).e)/2;

                linelist_up{i}(nline).len = sqrt( (linelist_up{i}(nline).s(1) - linelist_up{i}(nline).e(1))^2 + ...
                                            (linelist_up{i}(nline).s(2) - linelist_up{i}(nline).e(2))^2);    % link length


                linelist_up{i}(nline).ang = GetAngle((linelist_up{i}(nline).e(1) - linelist_up{i}(nline).s(1)), ...
                                    (linelist_up{i}(nline).e(2) - linelist_up{i}(nline).s(2))); % link orientation

                linelist_up{i}(nline).adj = [Inf Inf];

                linelist_up{i}(nline).nclq = 2; 

                linelist_up{i}(nline).prob = GetLineProb(ocsvm_clres_semi, linelist_up{i}(nline).s, linelist_up{i}(nline).e);

                linelist_up{i}(nline).nr = nline;
                
            end
        end
    end
    
    DrawLabelledLineList(ocsvm_clres_semi, linelist_up{i}, ones(size(linelist_up{i},2))); 
end
    
%%
for i=1:15
    if(~isempty(gt{i}))
        [linelist{i}, Cliques{i}, ~] = KovesiLineModel2MRFNetwork(rdclseglist{i}, clres{i}, 50, 0);
    end
end
%%
