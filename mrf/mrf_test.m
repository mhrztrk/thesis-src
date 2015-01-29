%% Preprocessing

% read image file
[imH imW imD] = size(img);

% import image from workspace

% create graph

% manually generate road training set
mask = zeros(imH, imW);  
while(1)
    mask_rd = roipoly(img(:,:,1:3));
    if(~isempty(mask_rd))
        mask = mask | mask_rd;
    end
    
    reply = input('do you want to continue?(Y/N): ', 's');
    
    if(strcmp(reply,'yes')||strcmp(reply,'Y')||strcmp(reply,'y'))
        break;
    end
end

% spectral transformation

img_med(:,:,1) = medfilt2(img(:,:,1));
img_med(:,:,2) = medfilt2(img(:,:,2));
img_med(:,:,3) = medfilt2(img(:,:,3));
img_med(:,:,4) = medfilt2(img(:,:,4));

img_ac_ao = imareaclose(imareaopen(uint8(img_med*255), 100, 8),100,8);

img_ac_ao_d(:,:,1) = double(img_ac_ao(:,:,1))/255;
img_ac_ao_d(:,:,2) = double(img_ac_ao(:,:,2))/255;
img_ac_ao_d(:,:,3) = double(img_ac_ao(:,:,3))/255;
img_ac_ao_d(:,:,4) = double(img_ac_ao(:,:,4))/255;

img_ac_ao_med(:,:,1) = medfilt2(img_ac_ao_d(:,:,1));
img_ac_ao_med(:,:,2) = medfilt2(img_ac_ao_d(:,:,2));
img_ac_ao_med(:,:,3) = medfilt2(img_ac_ao_d(:,:,3));
img_ac_ao_med(:,:,4) = medfilt2(img_ac_ao_d(:,:,4));

imgST = im_feature_transform(img_ac_ao_med, 'norm');

% generated road mask will be used to train road's spectral signature
clres = imclassify(imgST, mask, 'oc-svm');
figure;imshow(clres,[]);

% Apply Road template Matching filter to remove regions which have   
% structure not similar to roads. 
clresfilt = RoadTemplateMatchingFilterEx(clres, 5, 10, [7 16], 0);
figure;imagesc(clresfilt);

% skelotonization to obtain road centerlines.
th_bin = 0.16;
rdcline = bwmorph(clresfilt > th_bin ,'skel',Inf);
figure;imshow(rdcline,[]);

%% Kovesi - EdgeLinking Toolbox
    % Link edge pixels together into lists of sequential edge points, one
    % list for each edge contour. A contour/edgelist starts/stops at an 
    % ending or a junction with another contour/edgelist.
    % Here we discard contours less than 10 pixels long.
    [edgelist, labelededgeim] = edgelink(rdcline, 10);

    % Fit line segments to the edgelists
    tol = 2;  % Line segments are fitted with maximum deviation from
              % original edge of 2 pixels.

    % Road centerline segment list          
    rdclseglist = lineseg(edgelist, tol); 

    % Display the edgelists with random colours for each distinct edge 
    % in figure 2
    h = figure;
    imshow(img(:,:,1:3),[]);
    drawedgelist(rdclseglist, [imH imW], 1, 'rand', h); axis off 
    
    
    % deleting an element from a cell array
%     rdclseglist{12} = [];
%     rdclseglist(cellfun(@(rdclseglist) isempty(rdclseglist),rdclseglist))=[];
    
%%

% Find detected line segments
[linelist] = FindDetectedConnections(rdclseglist);

nlinedet = size(linelist,2);

% Find all possible candidate segments.
for i=1:2
    [linelist] = FindPossibleConnections(linelist);
    
    nlineconn(i) = size(linelist,2);

    [linelist] = ApplyLineCuttings(linelist);
    
    nlinecutt(i) = size(linelist,2);
end


%%
% second pass to remove duplicate lines

%%

% Create Adjacency Graph

nlink = size(linelist,2); % number of links
Adj = zeros(nlink,nlink);

for i=1:nlink
    for j=1:linelist(i).nconn
        Adj(i,linelist(i).conn(j)) = 1;
        Adj(linelist(i).conn(j),i) = 1;        
    end
end

% Determine the deflation angle between touching line pairs.  
for i=1:nlink
    for j=1:linelist(i).nconn
        linelist(i).defl(j) = defl_angle(linelist(i), linelist(linelist(i).conn(j)));
    end
end

% find maximal clique set
Cliques = FindCliques(linelist);

% Update line set structure "linelist"
[linelist] = UpdateLineInfo(linelist, Cliques, mat2gray(clres));

%%

figure; imshow(img(:,:,1:3),[]); hold on;
for i=1:nlinedet
    line([linelist(i).s(2) linelist(i).e(2)], [linelist(i).s(1) linelist(i).e(1)],...
        'Color',[0 0 1],'LineWidth',2);
    scatter(linelist(i).c(2), linelist(i).c(1), 'x', 'g');
end

for i=(nlinedet+1):nlineconn(1)
    line([linelist(i).s(2) linelist(i).e(2)], [linelist(i).s(1) linelist(i).e(1)],...
        'Color',[0 1 1],'LineWidth',2);
    scatter(linelist(i).c(2), linelist(i).c(1), 'x', 'g');
end

for i=(nlineconn(1)+1):nlinecutt(1)
    line([linelist(i).s(2) linelist(i).e(2)], [linelist(i).s(1) linelist(i).e(1)],...
        'Color',[1 0 0],'LineWidth',2);
    scatter(linelist(i).c(2), linelist(i).c(1), 'x', 'g');
end

for i=(nlinecutt(1)+1):nlineconn(2)
    line([linelist(i).s(2) linelist(i).e(2)], [linelist(i).s(1) linelist(i).e(1)],...
        'Color',[1 0 1],'LineWidth',2);
    scatter(linelist(i).c(2), linelist(i).c(1), 'x', 'g');
end

for i=(nlineconn(2)+1):nlinecutt(2)
    line([linelist(i).s(2) linelist(i).e(2)], [linelist(i).s(1) linelist(i).e(1)],...
        'Color',[1 1 0],'LineWidth',2);
    scatter(linelist(i).c(2), linelist(i).c(1), 'x', 'g');
end

%% Initalize Optimization parameters

nconn = size(Adj,1);    % # of connections

mrf.params.t1 = 0.2; 
mrf.params.t2 = 0.5;

mrf.params.Kl = 0.2; 
mrf.params.Ki = 0.2; 
mrf.params.Ke = 0.17; 
mrf.params.Kc = 0.2; 

mrf.params.eratio = (1/10);
mrf.params.allowjunc = 0;

labels = [ones(1,nlinedet) zeros(1,(nconn-nlinedet))]';

%% brute-force search for optimum solution

E_min = Inf;
labels_opt = labels;

for i=0:(2^nconn-1)    
    labels = de2bi(2^nconn-1- i, nconn)';
	E = CalculateEnergyLL(linelist, labels, Cliques, mrf);   
    if(E < E_min)
        E_min = E;
        labels_opt = labels;
        fprintf('Energy = %.3f\n', E_min);
    end
    
end


%% =~~=~~=~~=~~=~~=~~=~~=~~= Gibbs Sampler =~~=~~=~~=~~=~~=~~=~~=~~=~~=~~=~

step = 75;  % iteration count
T = 1;      % initial temperature

nLevel = 7; % # of temp levels (Simulated Annealing)

% initial graph labelling
%   simply assign 1 to detected line segmets, 0 to others
labels = [ones(1,nlinedet) zeros(1,(size(linelist,2)-nlinedet))]';

% MRF parameters

mrf = [];

mrf.sorm = 0;

mrf.params.Ke = 0.17; 
mrf.params.Kl = 0.13; 
mrf.params.Kc = 0.3; 
mrf.params.Ki = 0.2; 

mrf.params.allowjunc = 0;

[labels_GlbMin] = SimulatedAnnealingWithGibbsSampler(linelist, labels, step, nLevel, T, mrf, Cliques);

DrawLabelledLineList(img(:,:,1:3), linelist, labels_GlbMin);  

roadlines = GetFinalLineList(linelist, labels_GlbMin);

% final touch to connect the termination points
[LineList] = CompleteRoadNetwork(roadlines);
    

%% =~~=~~=~~=~~=~~=~~=~~=~~= END Gibbs Sampler =~~=~~=~~=~~=~~=~~=~~=~~=~~=          


%% =~~=~~=~~=~~=~~=~~=~~=~~=~~=~~= SORM =~~=~~=~~=~~=~~=~~=~~=~~=~~=~~=~~=~

[CC, CM, ST, pred] = sorm(im, 25, 15, 5);

%% =~~=~~=~~=~~=~~=~~=~~=~~= MRF Network (SORM) =~~=~~=~~=~~=~~=~~=~~=~~=~=

% Apply Road template Matching filter to remove regions which have   
% structure not similar to roads. 
clresfilt = RoadTemplateMatchingFilterEx(clres, 5, 10, [7 17], 0);
figure;imagesc(clresfilt);

th_bin = 0.17; 
rdcline = bwmorph(clresfilt > th_bin ,'skel',Inf);
figure;imshow(rdcline,[]);

[CC, CM, ST, pred] = sorm(rdcline, 25, 15, 5);

[linelist, Cliques] = SORMNetworkConstruct(CC, CM, clres_ac, 35, 0.2);

DrawLabelledLineList(img(:,:,1:3), linelist, ones(size(linelist,2))); 


%% =~~=~~=~~=~~=~~=~~=~~=~~= MRF Network (SORM) =~~=~~=~~=~~=~~=~~=~~=~~=~=

[linelist, Cliques] = SORMNetworkConstruct(CC, CM, mat2gray(clres > th_bin), 35, 0.3);

DrawLabelledLineList(img(:,:,1:3), linelist, ones(size(linelist,2))); 


%% =~~=~~=~~=~~=~~=~~=~~=~~= Gibbs Sampler(SORM) =~~=~~=~~=~~=~~=~~=~~=~~=~

nLevel = 7; % # of temp levels
step = 35;  % # of step in each level
c = 0.5;    % cooling schedule
Ti = 1;     % initial temperature

labels = ones(1,size(linelist,2))';

mrf = [];

mrf.sorm = 1;

mrf.params.Ke = 0.17; 
mrf.params.Kl = 0.13; 
mrf.params.Kc = 0.3; 
mrf.params.Ki = 0.2; 

mrf.params.eratio = (1/10);

mrf.params.allowjunc = 0;

[labels_GlbMin] = SAwithGibbsSampler(linelist, labels, step, nLevel, Ti, mrf, Cliques);

DrawLabelledLineList(img(:,:,1:3), linelist, labels_GlbMin);      

%%

Ef = CalculateEnergyLL(linelist, labels_GlbMin, Cliques, mrf);

extConn = [];

for i=1:size(Cliques,1)
    
    labels_tmp = labels_GlbMin;
    
    plist = Cliques{i}((labels_GlbMin(Cliques{i})==1));
    if(size(plist,2) == 1)
        nlist = Cliques{i}((labels_GlbMin(Cliques{i})==0));
        Er = zeros(size(nlist));
        for j=1:size(nlist,2)
            labels_tmp(nlist(j)) = 1;
            Er(j) = CalculateEnergyLL(linelist, labels_tmp, Cliques, mrf);
            labels_tmp(nlist(j)) = 0;
        end
    
        [val ind] = min(Er);
        if(val < Ef)
            extConn = [extConn nlist(ind)];
        end
    end
   
end

labelsExt = labels_GlbMin;
labelsExt(extConn) = ones(size(extConn)); 

%% =~~=~~=~~=~~=~~=~~=~~=~~= END Gibbs Sampler =~~=~~=~~=~~=~~=~~=~~=~~=~~= 


%% =~~=~~=~~=~~=~~=~~=~~= Optimization via ICM =~~=~~=~~=~~=~~=~~=~~=~~=~~=

step = 100;

E_new = Inf;
E_old = Inf;

ngrp = 4;
prob = [];

lprob = [linelist{5}.prob]; 
labels = (lprob > 0.5)';

nconn = size(labels,1);

for k=1:step
    for i=1:(nconn-(ngrp-1))
        labels_l = repmat(labels,1,2^ngrp);

        for m=1:2^ngrp
            labels_l(i:(i+ngrp-1),m)=de2bi((m-1),ngrp)';
        end            

       for m=1:2^ngrp
           prob(m) = 0;
           for n=1:ngrp
                prob(m) = prob(m) + exp(-MRFLocalEnergy(i+n-1, linelist{5}, labels_l(:,m), Cliques{5}, mrf));
           end
       end

        [~, ind] = max(prob);
        labels =  labels_l(:,ind);

        E = CalculateEnergyLL(linelist{5}, labels, Cliques{5}, mrf); 
        if(E < E_old)
            E_new = E; 
        end

        deltaE = abs(E_new - E_old);    % stop when energy change is small    
        E_old = E_new;

        fprintf('Energy = %.3f, deltaE = %.3f\n', E_new, deltaE);

    end
end
    
DrawLabelledLineList(img, linelist, labels);  

%% =~~=~~=~~=~~=~~=~~=~~=~~=~~= END ICM =~~=~~=~~=~~=~~=~~=~~=~~=~~=~~=~~=~

%%
figure; imshow(img,[]); hold on;

for i=1:size(linelist,2)
        line([linelist(i).s(2) linelist(i).e(2)], [linelist(i).s(1) linelist(i).e(1)],'Color',[0 0 1]);
        scatter(linelist(i).c(2), linelist(i).c(1), 'x', 'g');
end  
%%


