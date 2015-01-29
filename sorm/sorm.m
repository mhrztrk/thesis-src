%
%   Applies Self-Organizing Road Map Algorithm
%
%       I:  Binary Image, image pixels with value 1 indicate possible
%           road pixels.
%
%       d:  grid spacing
%
%       minDist        : minimum distance between the nodes 
%       minSampleCnt   : minimum number of samples within cluster 
%
function [CC, CM, ST, pred, Ic] = sorm(I, d, minDist, minSampleCnt, img)
    
    global dist;
    dist = d;

    % get image size
    [imH imW] = size(I);
    
    % find road pixels
    [X Y] = find(I == 1);
    m = zeros(size(X,1),2);
    m(:,1) = X;
    m(:,2) = Y;
    
    % initial cluster centers as uniform grid.  
    CCi = ClusterCentersInit(imH, imW, d); 
    
    % Apply k-medians algorithm %
    [CC,~,CM] = KMedians(m, CCi, minDist, minSampleCnt, I);
    
    Ic = ColorizeClusters([imH imW], CM);
    
    figure;imshow(Ic,[]); hold on;
    scatter(CC(:,2),CC(:,1),'r','o', 'filled', 'linewidth', 4);
%     
%     ST = 0;
%     pred = 0;
%     
%     VisualiseClusters(size(I), CC, CM);
%     
%     % Cluster Orientations
%     %  > co_pca : Cluster Orientations based on Principal Component Analysis
%     %  > co_ht  : Cluster Orientations based on Hough Transform
%     [co_pca co_ht] = ClusterOrintations(CM);
%     
%     % Deflection Angle Matrix 
%     D = GetDeflectionMatrix(CC, co_pca);
%     
%     % Proximity (Distance) Matrix
%     P = GetProximityMatrix(CC);
%     
%     % Read fuzzy inference system file (FIS)
%     fismat = readfis('LinkStrength.fis');
%     
%     % cluster count
%     cn = size(CC,1);
%     
%     DT = D';    
%     V1 = mat2gray(P(:));
%     V2 = mat2gray(D(:));
%     V3 = mat2gray(DT(:));
% 
%     % Determine link strengths using FIS
%     output = evalfis([V1 V2 V3],fismat);
%     Mfzzy = reshape(output,cn,cn);
%     Mfzzy = sparse(1 - Mfzzy);
%     
%     % Find minimum spanning tree based only on Proximity %
%     [ST,~] = graphminspantree(sparse(P));
%  
%     % Draw MST on original image %
%     %> ShowTreeOnImage(I, full(ST), CC);
%     
%     % Find minimum spanning tree based on Proximity and Deflection metrics %
%     [ST,pred] = graphminspantree(Mfzzy);
%     
%     % Draw MST on original image %
%     ShowTreeOnImage(img, full(ST), CC);
%     
%     %RNG = RNGraph(I, Mf, CC);
%     
%     %> ShowTreeOnImageSc(I, full(ST), Mfzzy, CC);
%     
%     scrsz = get(0,'ScreenSize');
%     set(gcf, 'Position',[1 1 scrsz(3) scrsz(4)]);
%     set(gcf, 'PaperPositionMode', 'auto');
%     saveas(gcf, sprintf('Minimum Spanning Tree (d=%2d,d_min=%2d).jpg', dist, minDist));
%         
end


%% 
% Calculates Deflection Angle Matrix 
% 
%   CC: Cluster Centers
%   CO: Cluster Orientation angles in degree
%
% deflection angle : Angle between link direction and principal 
% eigenvector direction
%
% D(i,j) -> Deflection angle between direction of link i to j and
%           principal eigenvector direction of cluster i.
%           D(i,j) and D(j,i) are not supposed to be equal.
%
function D = GetDeflectionMatrix(CC, CO)

    % cluster count
    cn = size(CC,1);

    D = zeros(cn,cn); 
    for i=1:cn
        for j=1:cn
            if(j ~= i)
                % link direction
                link_or = atand((CC(i,2)-CC(j,2))/(CC(i,1)-CC(j,1)));
                D(i,j) = abs(CO(i) - link_or);
            else
                D(i,j) = 90;
            end
        end
    end
end


%% 
% Calculates Proximity Matrix 
% 
%   CC: Cluster Centers
%
% Proximity(i,j) : euclidean distance between node i and j
%
%
function P = GetProximityMatrix(CC)

    P = pdist2(CC,CC);
    
end

%% 
function CC = ClusterCentersInit(imH, imW, d)
    
    crow_cnt = floor((imW - 1)/d);
    ccol_cnt = floor((imH - 1)/d);
    k = ccol_cnt*crow_cnt;
    CC = zeros(k, 2);
    
    CC(1:ccol_cnt,1) = d:d:(imH-1);
    CC(1:ccol_cnt,2) = d;
    
    for i=1:crow_cnt
        CC((ccol_cnt*(i-1)+1):(ccol_cnt*i),1) = CC(1:ccol_cnt,1);
        CC((ccol_cnt*(i-1)+1):(ccol_cnt*i),2) = d + (i-1)*d;
    end
    
end

%%
function VisualiseClusters(ISZ, CC, CM)

    Ic = ColorizeClusters(ISZ, CM);
    [co_pca co_ht] = ClusterOrintations(CM);
    
    figure;imshow(Ic,[]); hold on;
    scatter(CC(:,2),CC(:,1),'r','x');
    
    Vx = sind(co_ht); 
    Vy = cosd(co_ht);
    
    quiver (CC(:,2),CC(:,1),Vx',Vy','r','filled');
   
    figure;imshow(Ic,[]); hold on;
    scatter(CC(:,2),CC(:,1),'r','x');
    
    Vx = sind(co_pca); 
    Vy = cosd(co_pca);
    
    quiver (CC(:,2),CC(:,1),Vx',Vy','r','filled');
    
end

%%
function [co_pca co_ht] = ClusterOrintations(CM)
    
    % cluster count
    cn = size(CM,1);

    co_pca = size(cn,1);
    co_ht  = size(cn,1);
    
    for i=1:cn
        co_pca(i) = pca_orientation(CM{i});
        co_ht(i)  = ht_orientation(CM{i});
    end

end

%%
%   ISZ - Image Size [Height Width]
%   CM  - Cluster Members
%
function Ic = ColorizeClusters(ISZ, CM)
    
    % Pre-allocation
    Ic = zeros(ISZ(1), ISZ(2), 3);
    
    % cluster count
    cn = size(CM,1);
    
    % colored image size
    cim_sz = size(Ic);
    
    for i=1:cn
        
        % get points belongs to cluster i
        m = CM{i};
       
        % member count
        mcn = size(m,1);
        
        % assign a random color for each cluster
        Ic(sub2ind(cim_sz, m(:,1), m(:,2),   ones(mcn,1))) = rand(1)/2; 
        Ic(sub2ind(cim_sz, m(:,1), m(:,2), 2*ones(mcn,1))) = rand(1); 
        Ic(sub2ind(cim_sz, m(:,1), m(:,2), 3*ones(mcn,1))) = rand(1); 
        
    end
    
%     for i=1:im_sz(1)
%         
%         ind = find(g==i);
%         if ~isempty(ind)
%             % Generate a random color
%             re = rand(1); 
%             gr = rand(1); 
%             bl = rand(1); 
%             for j=1:im_sz(2)
%                 Ic(Y(ind(j)), X(ind(j)), 1) = re;
%                 Ic(Y(ind(j)), X(ind(j)), 2) = gr;
%                 Ic(Y(ind(j)), X(ind(j)), 3) = bl;
%             end
%         end
%     end

end


%% 
function [Vp Ir R] = ClusterEig(I, X, Y, c, g)
    %
    % |A  B|
    % |B  C|
    %
    % A = x^2, B = xy, C = y^2
    %
    R = zeros(size(c,1),1);
    Ir = zeros(size(I));
    Vp = zeros(size(c,1),2);
    
    for i=1:size(c,1)
        ind = find(g==i);
        if ~isempty(ind)
            A = 0;
            B = 0;
            C = 0;
            for j=1:size(ind,1)
                A = A + (X(ind(j))-c(i,2))^2; 
                B = B + (X(ind(j))-c(i,2))*(Y(ind(j))-c(i,1));
                C = C + (Y(ind(j))-c(i,1))^2;
            end
            
            [V D] = eig([A B;B C]);
            lambda = abs(diag(D));
            
            [mval ind] = max(lambda);
            R(i) = mval / min(lambda);
            Vp(i,:) = V(:,ind)'; 
            
            for j=1:size(ind,1)
                Ir(Y(ind(j)), X(ind(j))) = R(i);
            end
            
        end
    end
end

%%
function d=DistMatrix(A,B)
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % DISTMATRIX return distance matrix between points in A=[x1 y1 ... w1] and in B=[x2 y2 ... w2]
    % Copyright (c) 2005 by Kardi Teknomo,  http://people.revoledu.com/kardi/
    %
    % Numbers of rows (represent points) in A and B are not necessarily the same.
    % It can be use for distance-in-a-slice (Spacing) or distance-between-slice (Headway),
    %
    % A and B must contain the same number of columns (represent variables of n dimensions),
    % first column is the X coordinates, second column is the Y coordinates, and so on.
    % The distance matrix is distance between points in A as rows
    % and points in B as columns.
    % example: Spacing= dist(A,A)
    % Headway = dist(A,B), with hA ~= hB or hA=hB
    %          A=[1 2 3; 4 5 6; 2 4 6; 1 2 3]; B=[4 5 1; 6 2 0]
    %          dist(A,B)= [ 4.69   5.83;
    %                       5.00   7.00;
    %                       5.48   7.48;
    %                       4.69   5.83]
    %
    %          dist(B,A)= [ 4.69   5.00     5.48    4.69;
    %                       5.83   7.00     7.48    5.83]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    [hA,wA]=size(A);
    [hB,wB]=size(B);
    if wA ~= wB,  
        error(' second dimension of A and B must be the same'); 
    end
    
    C = zeros(1,wa);
    D = zeros(1,wa);
    
    for k=1:wA
      C{k}= repmat(A(:,k),1,hB);
      D{k}= repmat(B(:,k),1,hA);
    end
    S=zeros(hA,hB);
    for k=1:wA
      S=S+(C{k}-D{k}').^2;
    end
    d=sqrt(S);
end

%%
function ind = GetMemberInd2(A, B, blocksz)
    
    [hA,wA]=size(A);
    [hB,wB]=size(B);
    if wA ~= wB,  
        error(' second dimension of A and B must be the same'); 
    end
    
    C = cell(1,wA);
    D = cell(1,wB);
    c = idivide(int32(hA),blocksz);
    rm = mod(hA, blocksz);
    ind = zeros(hA,1);
    
    for i=1:c

        s_ind = (blocksz*(i-1)+1);
        e_ind = blocksz*i;
        
        AA = A(s_ind:e_ind,:); 
        
        for k=1:wA  
          C{k} = repmat(AA(:,k), 1, hB);
          D{k} = repmat( B(:,k), 1, blocksz);
        end
        
        S=zeros(blocksz,hB);

        for k=1:wA
          S=S+(C{k}-D{k}').^2;
        end
        
        [~, ind(s_ind:e_ind,1)] = min(S, [], 2);
        
    end
    
    if (rm ~= 0)
        s_ind = (blocksz*(c-1)+1);
        e_ind = rm;
        blocksz_rm = e_ind-s_ind+1; 
        
        AA = A(s_ind:e_ind,:); 
        
        for k=1:wA  
          C{k} = repmat(AA(:,k), 1, hB);
          D{k} = repmat( B(:,k), 1, blocksz_rm);
        end
        
        S=zeros(blocksz_rm, hB);

        for k=1:wA
          S=S+(C{k}-D{k}').^2;
        end
        
        [~, ind(s_ind:e_ind,1)] = min(S, [], 1);
        
    end
    
end

%%
function ind=GetMemberInd(A,B)

    [hA,wA]=size(A); 
    [~,wB]=size(B);
    if wA ~= wB,  
        error(' second dimension of A and B must be the same'); 
    end
    
    ind = zeros(hA,1);
  
    for i=1:hA
        D = (B(:,1) - A(i,1)).^2 + (B(:,2) - A(i,2)).^2;
        [~, ind(i,1)] = min(D);
    end
    
end

%%
% CC - Cluster Centers
% ID - Cluster ID of each point
% CM - Members of each point
%
function [CC ID CM]= KMedians(m, sp, min_dist, minSampleCnt, I)

    temp=zeros(size(m, 1),1);   % initialize as zero vector
    k = size(sp, 1);
    c = sp;
    iter = 0;
    
    global dist;
    
    % Show initial centroid locations
    figure;imshow(I);hold on
    scatter(sp(:,2),sp(:,1),'b','o', 'filled');     
    
    figure;imshow(I);hold on
    voronoi(sp(:,2),sp(:,1), ':');
    scatter(sp(:,2),sp(:,1),'b','o', 'filled');     
    
    while 1
        ID = GetMemberInd(m, c);  % calculate objcets-centroid distances
        
        if ID == temp
            % remove centroids with no member
            j=0;
            c2 = c;
            ds = zeros(k,1);
            for i=1:k
                f = find(ID==i);
                if f          
                    ds(i,1) = 1;
                    j = j + 1;
                    c2(j,:) = c(i,:);
                else
                    ds(i,1) = 0;
                end
            end
            
            if (j ~= 0)
                c = c2(1:j,:);
            else
                c = 0;
            end
            
            % Call GetMemberInd one more time to update g, 
            ID = GetMemberInd(m, c);
            
            break;          % stop the iteration
        else
            temp = ID;         % copy group matrix to temporary variable
        end
        
        
        %% Remove clusters which have samples less than min_sample_th
        
%         c_up = c;
%         j = 0;
%         
%         for i=1:size(c,1)
%             f = find(ID==i);
%             if size(f,1) >= min_sample_th          
%                 j = j + 1;
%                 c_up(j,:) = c(i,:);
%             end
%         end
%         c = c_up(1:j,:);
        
        %% Draw and Save
%         c_w_member = zeros(size(c));
%         c_wo_member = zeros(size(c));
%         j = 0;
%         k = 0;
%         for i=1:size(c,1)
%                 f = find(ID==i);
%                 if size(f,1) > 0          
%                     j = j + 1;
%                     c_w_member(j,:) = c(i,:);
%                 else
%                     k = k + 1;
%                     c_wo_member(k,:) = c(i,:);
%                 end
%         end
%         c_w_member = c_w_member(1:j,:);
%         c_wo_member = c_wo_member(1:k,:);
%     
%         iter = iter + 1;
%         
%         % Draw result
%         figure;imshow(I);hold on
%         scatter(c_w_member(:,2),c_w_member(:,1),'r','o', 'filled'); 
%         scatter(c_wo_member(:,2),c_wo_member(:,1),'b','o', 'filled');  
%         
%         % Save result
%         scrsz = get(0,'ScreenSize');
%         set(gcf, 'Position',[1 1 scrsz(3) scrsz(4)]);
%         set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf, sprintf('KMedians(d=%2d,d_min=%2d)_iter%02d.jpg', dist, min_dist, iter));
%         close gcf;  
        
        %%
        
        for i=1:k
            f = find(ID==i);
            
            % If number of samples belongs to node i is below the specified
            % threshold min_sample_the, delete the node.
            if f            % only compute centroid if f is not empty
                c(i,:) = median(m(ID==i,:),1);
                % c(i,:) = mean(m(ID==i,:),1);
            end
            
        end

        % Node Pruning: Merge the node pairs if the ecludian distance 
        % between them is below the min_dist.
        % 
        D = GetDistMatrix(c);

        % Check if there are nodes to be merged.
        if (~isempty(find(D < min_dist, 1)))

            D(D >= min_dist) = Inf;
            c_pruned = c;

            while(1)
                % find the node pair in which the distance between the nodes is minimum. 
                [v, ind] = min(D(:));
                if (~isinf(v))
                    [~, jj] = ind2sub(size(D), ind);

                    % remove one of the node in the pair by replacing the location of it by 0 
                    c_pruned(jj,:) = [0 0];

                    % update the Distance matrix
                    % The Node jj is removed, so delete the entries related
                    % with node jj.
                    D(:,jj) = Inf;
                    D(jj,:) = Inf;
                else
                    c = c_pruned(c_pruned(:,1)~=0 & c_pruned(:,1)~=0,:);
                    break
                end

            end

        end
        
                   
    end
    
    % cluster count
    cn = size(c,1);
    
    % Construct Cluster structure
    
    k = 0;
    for i=1:cn 
        S = m(ID==i,:);
        if(size(S,1) >= minSampleCnt)
            k = k + 1;
            CM{k} = S;
            CC(k,:) = round(c(i,:)); 
        end
    end
    CM = CM';
    
    % Do we really need to round cluster centers???????????????????????????
    
   figure;imshow(I,[]); hold on;
   scatter(CC(:,2),CC(:,1),'b','o', 'filled');
    
    figure;imshow(I);hold on
    voronoi(CC(:,2),CC(:,1),':');
    scatter(CC(:,2),CC(:,1),'b','o', 'filled');  
   
end

%%
function [D] = GetDistMatrix(c)

    cn = size(c, 1);
    
%     X1 = repmat(c(:,1),cn,1);
%     X2 = repmat(c(:,1),1,cn);
%     X2 = X2';
%     X2 = X2(:);
%     
%     Y1 = repmat(c(:,2),cn,1);
%     Y2 = repmat(c(:,2),1,cn);
%     Y2 = Y2';
%     Y2 = Y2(:);
    
% Improvement
    X2 = repmat(c(:,1),1,cn);
    X1 = X2(:);
    X2 = X2';
    X2 = X2(:);
    
    Y2 = repmat(c(:,2),1,cn);
    Y1 = Y2(:);
    Y2 = Y2';
    Y2 = Y2(:);


    D = sqrt((X1 - X2).^2 + (Y1 - Y2).^2);
    D = reshape(D, cn, cn);
    
    % we need only lower triangular matrix 
    D = tril(D);
    
    % it is required to initialize diagonal and upper triangular part of 
    % the D matrix with Inf. 
    D(D==0) = Inf;

% Older Method to calculate Distance Matrix, 
% the method above is much more faster than this.   
%
% D = Inf(cn,cn);    
% for i=2:cn
%     for j=1:(i-1)
%         D(i,j) = sqrt((c(i,1) - c(j,1))^2 + (c(i,2) - c(j,2))^2);
%     end
% end

end



%%
function cnt = RoadPixelCnt(I, c1, c2)
        [Y,X] = GetPixelsOnTheLine(c1,c2);
        cnt = 0;
        n = size(X,2);
        for i=1:n
            cnt = cnt + I(X(1,i),Y(1,i));
        end
        %msk = (I(X,Y) == 1);
        %cnt = sum(msk(:));
end

%%
function [X,Y] = GetPixelsOnTheLine(c1,c2)

    % BRESENHAM: Generate a line profile of a 2d image 
    %            using Bresenham's algorithm
    % [myline,mycoords] = bresenham(mymat,mycoordinates,dispFlag)
    %
    % - For a demo purpose, try >> bresenham();
    %
    % - mymat is an input image matrix.
    %
    % - mycoordinates is coordinate of the form: [x1, y1; x2, y2]
    %   which can be obtained from ginput function
    %
    % - myline is the output line
    %
    % - mycoords is the same as mycoordinates if provided. 
    %            if not it will be the output from ginput() 
    % Author: N. Chattrapiban
    %
    % Ref: nprotech: Chackrit Sangkaew; Citec
    % Ref: http://en.wikipedia.org/wiki/Bresenham's_line_algorithm
    % 
    % See also: tut_line_algorithm
    
    x(1) = round(c1(1,1));
    x(2) = round(c2(1,1));
    y(1) = round(c1(1,2));
    y(2) = round(c2(1,2));
    
    steep = (abs(y(2)-y(1)) > abs(x(2)-x(1)));

    if steep 
        [x,y] = swap(x,y); 
    end

    if x(1)>x(2), 
        [x(1),x(2)] = swap(x(1),x(2));
        [y(1),y(2)] = swap(y(1),y(2));
    end

    delx = x(2)-x(1);
    dely = abs(y(2)-y(1));
    error = 0;
    x_n = x(1);
    y_n = y(1);
    
    if y(1) < y(2) 
        ystep = 1; 
    else
        ystep = -1;
    end
    
    X = zeros(1,delx+1);
    Y = zeros(1,delx+1);
    
    for n = 1:(delx+1)
        if steep
            X(n) = x_n;
            Y(n) = y_n;
        else
            X(n) = y_n;
            Y(n) = x_n;
        end    
        x_n = x_n + 1;
        error = error + dely;
        if bitshift(error,1) >= delx % same as -> if 2*error >= delx, 
            y_n = y_n + ystep;
            error = error - delx;
        end    
    end

end

%%
function [q,r] = swap(s,t)
    % function SWAP
    q = t; r = s;
end

%%
function ShowTreeOnImage(I, M, CC)

    % display image 
    figure, imshow(I, []), hold on;
   
    % display minumum spanning tree
    w = size(M,1);
    for i=1:w
        ind = find(M(:,i) ~= 0);
        if ~isempty(ind)
            for j=1:size(ind,1)
                line([CC(i,2) CC(ind(j),2)],[CC(i,1) CC(ind(j),1)],...
                    'Color',[0 0 1],'LineWidth',4);              
            end
        end
    end
  
    % display cluster center locations found by K-medians algorithm
   scatter(CC(:,2),CC(:,1), 'r','o','filled');  
    
end


%%
function ShowTreeOnImageSc(I, M, LS, CC)

    % display image 
    figure, imshow(I, []), hold on;
    
    % display cluster centroids which are found previously by applying 
    % K-medians algorithm. 
    scatter(CC(:,2),CC(:,1),'r','x','LineWidth', 2);
    
    ls_min = min(M(M ~= 0));
    ls_max = max(M(M ~= 0));
    ls_dst = (ls_max - ls_min) / 2;
    ls_mid = ls_min + ls_dst;
    
    
    % display minumum spanning tree
    w = size(M,1);
    for i=1:w
        ind = find(M(:,i) ~= 0);
        if ~isempty(ind)
            for j=1:size(ind,1)
                
                val = M(ind(j), i);
                
                if (val > ls_mid)
                    b = (val - ls_mid) / ls_dst;
                else
                    b = 0;
                end
                
                if (val < ls_mid)
                    g = (val - ls_min) / ls_dst;
                else
                    g = (ls_max - val) / ls_dst;
                end
                
                if (val < ls_mid)
                    r = (ls_mid - val) / ls_dst;
                else
                    r = 0;
                end
                
                line([CC(i,2) CC(ind(j),2)],[CC(i,1) CC(ind(j),1)],'Color',[r g b],'LineWidth',4);
            end
        end
    end
end

%%
function RNG = RNGraph(I, M, c)

    % number of nodes
    n = size(M,1);
    
    % number of node pairs
    np = n*(n-1)/2;
    
    % Make diagonal entries Inf
    M = M + diag(inf(1,size(M,1)));
    
    RNG = zeros(size(M));
    
    for i=1:(n-2)
        for j=(i+1):n
            
            not_rng = 0;
            
            for k=1:n 
                if (k~=i && k~=j)
                    if(M(i,j) >= M(i,k)) && (M(i,j) >= M(j,k))
                        not_rng = 1;
                        break;
                    end
                end
            end
            
            if (not_rng == 0)
                RNG(i,j) = 1;
            end
            
        end
    end

    [Y X] = find(RNG == 1);
    
    % display image 
    figure, imshow(I, []), hold on;
    
    % display cluster center locations which are found previously by applying 
    % K-medians algorithm. 
    scatter(c(:,2),c(:,1),'r','x');
    
    % display minumum spanning tree
    for i=1:size(Y,1)
        line([c(Y(i),2) c(X(i),2)],[c(Y(i),1) c(X(i),1)]);
    end
    
end

%
%   == Principal Component Analysis Orientation ==
% 
% V -> Voronoi set: image pixels associated with node n 
% n -> node 
%                         2 * sigma_xy  
% ang = 0.5 * arctan( --------------------- );
%                     sigma_x^2 - sigma_y^2
%
%   sigma_x, sigma_y -> standard deviations in the direction of x and y
%   sigma_xy         -> covariance
%
%   |--------> Y 
%   |  
%   |  Image
%   |
%  \./
%   X

function ang = pca_orientation(V, n)
    
    % samples 
    X = V(:,1);
    Y = V(:,2);
    N = size(V,1);
    
    % sample mean
    m_x = sum(X(:))/N;
    m_y = sum(Y(:))/N;

    % sample variance
    var_x = ((X-m_x)'*(X-m_x))/(N*N);
    var_y = ((Y-m_y)'*(Y-m_y))/(N*N);
    
    % covariance
    cov_xy = ((X-m_x)'*(Y-m_y))/(N*N);
    
    %
    %     | sigma_x^2   sigma_xy |
    % C = |                      |
    %     | sigma_xy   sigma_y^2 |
    %
    C = [var_x cov_xy; cov_xy var_y];
    
    % determine eigenvectors and eigenvalues
    [VV DD] = eig(C);
    
    % orientation of cluster determined by the eigenvector corresponding to
    % maximum eigenvalue.
    [~, mx_eig_ind] = max(diag(DD)); 
    
    % calculate orientation angle in degree
    ang = atand(VV(mx_eig_ind,2)/VV(mx_eig_ind,1));
    
end


%
%   == Hough Transform Orientation ==
%
% V -> Voronoi set: image pixels associated with node n 
% n -> node 
%
function ang = ht_orientation(V, n)
    
    [h ,~] = max(V(:,1));
    [w, ~] = max(V(:,2));
    
    % reconstruct image from Voronoi set
    ind = sub2ind([h w], V(:,1), V(:,2));
    bw = zeros(h,w);
    bw(ind) = 1;
    
    % determine peak point of Hough Transform
    [H,~,~] = hough(bw);
    P  = houghpeaks(H);
    ang = P(2);
    
end


