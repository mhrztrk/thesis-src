%%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
%   Thesis Work:
%
%       Markov Random Field based Road Network Extraction from High Resolution
%       Satellite Images 
%   
%   Mahir OZTURK / mhrztrk@gmail.com
%
%%=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
%% Algorithm Parameters

% params.prefilt = 1;                 % include pre-filtering stage 
params.prefilt = 0;               % do not include pre-filtering

params.seedsel = 'auto';            % automatic seed point selection via ACE
    %params.autoseeduse = 'single';    % Use centerline segment seperately, then fuse classifications
    params.autoseeduse = 'all';     % Use all segments for classification
    % params.dilateseedpoints = 1;    % Thicken apar edges  
    
% params.seedsel = 'manual';        % manual seed point selection

params.linemodel = 'sorm';          % use SORM algorithm to generate candidate road segments    
% params.linemodel = 'skel';        % use skeloton extraction

params.clsmethod = 'ocsvm';         % classification method selection 
% params.clsmethod = 'gmm';

% MRF Relation method
params.mrfrelax = 'gibbs';          % Simulated Anneaking with Gibbs sampler
% params.mrfrelax = 'icm';          % Iterated conditional mode

% specify an image name
img_name = 'img';

% directory for saving results
base_dir = '';


%% Function search directories


%% Pre-filtering (Optional) / Noise (or detail) removing

if(params.prefilt)
    %  Bilateral Filter (to obtain cartoonish images)

    % # of spectral bands (or features) in image 
    nband = size(img,3);
    img_bf = zeros(size(img)); 

    if(nband == 3)
        % Apply color bilateral filtering.
        img_bf = bfilter2(img, 5, [3 0.05]);
    else
        for j=1:nband
            % Apply grey-scale bilateral filtering on each band
            img_bf(:,:,j) = bfilter2(img(:,:,j), 5, [3 0.05]);
        end
    end
end

%% Spectral Transformation / Feature Extraction

% # of spectral bands in image
nband = size(img,3);

% Principal Component Transform
if(params.prefilt)
    model_pca = imgpca(uint8(mat2gray(img_bf)*255),... % scale image
                            'auto', nband);
else
    model_pca = imgpca(uint8(mat2gray(img)*255),... % scale image
                            'auto', nband);
end

% Take first two components as new feature set (first two spectral is generally 
% enough for representing whole spectral content of an 3 or 4 band images.) 
% Scale feature set between -1 and 1 for classification
fset(:,:,1) = mat2gray(model_pca.scores(:,:,1))*2 - 1;  % PC1
fset(:,:,2) = mat2gray(model_pca.scores(:,:,2))*2 - 1;  % PC2

%% Road Seed point selection

[imH imW imD] = size(img);


if(imD >= 3)
    img_disp = img(:,:,1:3);            % image for displaying results
    img_gray = rgb2gray(img(:,:,1:3));  % gray-level image for edge detection
else
    % if # of spectral bands is less than three use first band for edge
    % detection.
    img_disp = img(:,:,1);              % image for displaying results
    img_gray = img(:,:,1);              % gray-level image for edge detection
end


if(strcmp(params.seedsel,'auto'))
    % automatic seed selection (ACE)
    img_ace = ace( img_gray, ... % gray-level or pan image
                        5,  ...  % minimum distance between apar-edge points
                        25, ...  % maximum distance 
                        10, ...  % maximum angle between gradient orientations    
                        0, ...  % minimum apar-edge length   
                        0,  ...  % save intermadiate results (see ace.m)   
                        img_name);
    
    % thicken apar-edges to increase training set.
    % TODO: Use the region between the road edges for training.
    if(params.dilateseedpoints)
        img_ace = bwmorph(img_ace, 'dilate');
    end

    img_ace = impyramid(img_ace, 'expand');
    
    ace_flt = RoadTemplateMatchingFilterEx2( img_ace, 3, 10, [3 5], 0, 3);
    img_ace = cc_threshold(ace_flt > 0.3,100,0);
    
    img_ace = impyramid(img_ace, 'reduce');
    mask = img_ace;
    
    figure;imshow(img_ace);
else
    % manual seed selection
    mask = zeros(imH, imW);  
    while(1)
        mask_rd = roipoly(img(:,:,1:3));
        if(~isempty(mask_rd))
            mask = mask | mask_rd;
        end

        reply = input('Press Q to exit: ', 's');

        if(strcmp(reply,'Q')||strcmp(reply,'q'))
            break;
        end
    end
end


%% Classification (Road model generation)
clres = [];
qmetrics = [];

if(strcmp(params.clsmethod,'ocsvm'))
    nu = 0.01;
else
    nu = [];
end
range = [-1 1]; 

if(strcmp(params.seedsel,'auto') && strcmp(params.autoseeduse,'single'))
    
        % GEnerate a seperate road model for each apar edge
        CC = bwconncomp(img_ace);
        
        for i=1:CC.NumObjects
            mask = zeros(imH, imW);
            mask(CC.PixelIdxList{i}) = 1;
            [clres{i} qmetrics{i}] = imclassify( fset, ...
                                mask, ...
                                zeros(imH, imW), ...            % ground truth for performance evaluation
                                params.clsmethod, ...           % training method
                                nu, ...                         % nu-parameter of OCSVM classifier
                                range, ...                  
                                img_name, ...        
                                sprintf('%s/clres/%d', base_dir,i)); % directory for saving results
        end
else
    
    % Generate a unique road model from all apar edges 
    [clres qmetrics] = imclassify( fset, ...
                        img_ace, ...
                        zeros(imH, imW), ...            % ground truth for performance evaluation
                        params.clsmethod, ...           % training method
                        nu, ...                         % nu-parameter of OCSVM classifier
                        range, ...                  
                        img_name, ...        
                        sprintf('%s/clres/%d', base_dir,i)); % directory for saving results         
end

%% Road Class Refinement

% initialization
clres_filt = [];
clres_filt_th = zeros(imH, imW);

% threshold depends on iteration cound and min/max road width
filt_th = 0.11; 

% minimum # of pixels threshold in connected component
min_cc = 100;

if(strcmp(params.seedsel,'auto') && strcmp(params.autoseeduse,'single'))
    
    for i=1:size(clres,2)
        clres_filt{i} = RoadTemplateMatchingFilterEx2( clres{i},  ...     % binary or gray-level image 
                                                    5,      ...     % itertion count
                                                    10,     ...     % angle between rotations
                                                    [7 13], ...     % min&max road width    
                                                    0,      ...     % display intermadiate results
                                                    3);             % method
        % convert to binary image and fuse results
        clres_filt_th = clres_filt_th | cc_threshold(clres_filt{i} > filt_th, min_cc, 0);
    end
    
    figure; imshow(clres_filt_th, []);
        
else
    clres_filt = RoadTemplateMatchingFilterEx2( clres,  ...     % binary or gray-level image 
                                                5,      ...     % itertion count
                                                10,     ...     % angle between rotations
                                                [7 13], ...     % min&max road width    
                                                0,      ...     % display intermadiate results
                                                3);             % method
    % convert to binary image
    clres_filt_th = (clres_filt > filt_th);
    
    figure;
    subplot(1,2,1),imshow(clres_filt, []);
    subplot(1,2,2),imshow(clres_filt_th, []);
end

%% Road Centerline Extraction

% skeloton extraction
img_skel = bwmorph(clres_filt_th,'skel',Inf);

% prune branches with lenght less than or equal to 15
img_bp = bwmorph(img_skel,'branchpoints');
img_bp = bwmorph(img_bp,'dilate');
img_skel(img_bp==1) = 0;
img_skel = cc_threshold(img_skel, 15, 0);

figure;imshow(bwmorph(img_skel,'dilate'),[]);

clear img_bp

if (strcmp(params.linemodel, 'sorm'))
    
    clres_filt_th_ref = bwmorph(img_skel,'dilate');
    
    [CC, ~ ] = sorm(clres_filt_th_ref, ...     
                                      25, ...     % grid spacing
                                      12, ...     % min distance between cluster centroids
                                      5);         % min number of samples in a cluster (i.e. clusters with samples less than this will be removed) 
                                        
                                        
   [linelist, Cliques] = SORMNetworkConstruct( CC,         ... % Cluster centroids
                                               [],         ... 
                                               clres,   ... % classified image
                                               50,         ... % max length of connecting line segment 
                                               0.1);           % min line probability 
     
else
    % Kovesi - EdgeLinking Toolbox
    [edgelist, labelededgeim] = edgelink(img_skel, 10);

    % Fit line segments to the edgelists
    tol = 2;  % Line segments are fitted with maximum deviation from
              % original edge of 2 pixels.

    % Road centerline segment list          
    rdclseglist = lineseg(edgelist, tol); 

    % Display the edgelists with random colours for each distinct edge 
    % in figure 2
    h = figure;
    imshow(img(:,:,1:3),[]);
    drawedgelist(rdclseglist, [size(img,1) size(img,2)], 1, 'rand', h); axis off  
    
    [linelist, Cliques, ~] = KovesiLineModel2MRFNetwork(  rdclseglist, ... % road segment list obtained by lineseg() func
                                                                clres,       ... % classified image  
                                                                50,          ... % max length of connecting line segment   
                                                                0.1);            % min line probability   
    
end


%% Road Network Formation

nLevel = 7; % # of temp levels
step = 50;  % # of step in each level
c = 0.5;    % cooling schedule
Ti = 2;     % initial temperature

mrf = [];
mrf.sorm = 1;

mrf.params.Ke = 0.3;   
mrf.params.Kl = 0.01; % Since lengths of segments are similar, Kl parameter is kept so small.  
mrf.params.Kc = 0.3; 
mrf.params.Ki = 0.4; 

mrf.params.Vo = 0.7;

mrf.params.eratio = (1/25);
mrf.params.allowjunc = 0;

% save connection info for debugging
save_linelist_info(linelist, sprintf('linelist_info_%s.txt', img_name));  
    
if (strcmp(params.mrfrelax,'gibbs'))    
    % MRF Relaxation / Simulated Annealing with Gibbs Sampler %% 
    
    % initial labeling / not important for stochastic relaxation
    iLabels = ones(1,size(linelist,2))';

    % optimization
    opt_labels = SAwithGibbsSampler(linelist, iLabels, step, nLevel, Ti, mrf, Cliques, 0);
    
else
    % MRF Relaxation / ICM %%
    
    % initial labeling
    iLabels = ([linelist.prob] > 0.5)';

    opt_labels = MRFRelaxWithICM(linelist, iLabels, step, mrf, Cliques);
end

% Dispay result 
DrawLabelledLineList(img_disp, linelist, opt_labels);

%% Performance evaluation
    
if(~exist('img_gt'))

    ext = ones(imH,imW);

    % generate extracted road network
    for j=1:size(linelist,2)
        if(opt_labels(j) == 1)
            [~,~,ext] = bresenham(ext,[ linelist(j).s(2), ...
                                           linelist(j).s(1); ...
                                           linelist(j).e(2),...
                                           linelist(j).e(1)],0);          
        end
    end
    ext = imcomplement(ext);   
    figure;imshow(ext,[]);
    
    % reference road network
    ref = bwmorph(gt,'skel',Inf);
    
    % thickened extracted & reference road networks
    ref_buf = bwmorph(ref,'dilate',5);
    ext_buf = bwmorph(ext,'dilate',5);

    % matched reference network
    match_ref = ext & ref_buf;

    % matched extracted network
    match_ext = ref & ext_buf;
    
    % calculate quality metrics
    completeness = sum(match_ref) / sum(ref);
    correctness  = sum(match_ext) / sum(ext);
    quality      = completeness * correctness / (completeness + completeness * correctness - correctness);
   
    % visualize evaluation result
    eval = zeros(imH,imW,3);
    % missed extractions
    R = bwmorph(ref & imcomplement(ext_buf), 'dilate');
    % correct extractions
    G = bwmorph(match_ext,'dilate');
    % false extractions
    B = bwmorph(ext & imcomplement(ref_buf),'dilate');
    % empty locations
    W = imcomplement(R) & imcomplement(G) & imcomplement(B);

    eval(:,:,1) = R + 0.5 * W;
    eval(:,:,2) = G + 0.5 * W;
    eval(:,:,3) = B + 0.5 * W;

    figure;imshow(eval,[])

    imwrite(eval, sprintf('%s/rnet/%s_eval_result.png',base_dir,img_name));
else
    msgbox('Ground Truth is not found, evaluation will not be performed!', 'Warning');
end
%%


