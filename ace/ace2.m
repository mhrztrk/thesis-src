function [ Iout ] = ace2( img, min_dist, max_dist, min_defl, cc_th, display, img_name)
%
% Anti-parallel-edge Centerline Extraction (ACE) Algorithm
%

    img = mat2gray(img);

    % Image size
    [imH imW] = size(img);
    
    % find canny edges
    % BW = edge(I,'canny',thresh,sigma)
    CannyEdge = edge(img, 'canny');

    % 2-D 3x3 Sobel Masks
    sb_x = fspecial('sobel');
    
%     sb_x = [1   2  0  -2   1 ;...
%             4   8  0  -8  -4 ;...
%             6  12  0 -12  -6 ;...
%             4   8  0  -8  -4 ;...
%             1   2  0  -2  -1];
        
    sb_y = imrotate(sb_x, 90, 'crop');
        
    % Apply Sobel Filter to determine the edge orientation
    grad_x= imfilter(img, sb_x);
    grad_y= imfilter(img, sb_y);
    
    CL.ace{1} = zeros(imH, imW); 
    CL.ace{2} = zeros(imW, imH);
  
    CL.CannyEdge{1} = CannyEdge;
    CL.CannyEdge{2} = CannyEdge';
    
    CL.GradX{1} = grad_x;
    CL.GradX{2} = grad_x';
    
    CL.GradY{1} = grad_y;
    CL.GradY{2} = grad_y';
   
    CL.GradAngle{1} = GetAngleM(grad_x, grad_y);
    CL.GradAngle{2} = CL.GradAngle{1}'; 
     
    CL.AparEdgeLocX{1} = zeros(imH, imW,2);
    CL.AparEdgeLocY{1} = zeros(imH, imW,2);
    
    CL.AparEdgeLocX{2} = zeros(imH, imW,2);
    CL.AparEdgeLocY{2} = zeros(imH, imW,2);
    
    max_dist = max_dist * 2;
    
    % search centerline points in x and y direction
    for d=1:2
        
        % d = 1 -> search in y direction
        % d = 2 -> search in x direction
        
        for i=1:size(CL.CannyEdge{d},2)

            ce_line = CL.CannyEdge{d}(:,i);
            ep = find(ce_line == 1);
            LP = size(CL.CannyEdge{d},1);
            
            if(isempty(ep))
                continue;
            end

            for j=1:size(ep,1)

                ind = ep(j);
                ce_line = CL.CannyEdge{d}(:,i);
                
                % Apply dmin and dmx thresholds.
                if ~(ind <= max_dist)
                    ce_line(1:(ind-max_dist)) = 0;
                end

                if ~((ind+max_dist) > LP)
                    ce_line((ind+max_dist):end) = 0;
                end

                if ~(ind <= min_dist)
                    if ~((ind+min_dist) > LP)
                        ce_line((ind-min_dist):(ind+min_dist)) = 0;
                    else
                        ce_line((ind-min_dist):end) = 0;
                    end
                else
                    ce_line(1:(ind+min_dist)) = 0;
                end

                pep = find(ce_line == 1);

                if(isempty(pep))
                    continue;
                end

                %ang_j = GetAngle(CL.GradX{d}(ind,i), CL.GradY{d}(ind,i));
                ang_j = CL.GradAngle{d}(ind,i);
                
                for k=1:size(pep,1)

                    %ang_k = GetAngle(CL.GradX{d}(pep(k), i), CL.GradY{d}(pep(k), i));
                    ang_k = CL.GradAngle{d}(pep(k), i);
                    
                    if ( abs(abs(ang_k - ang_j) - 180) < min_defl) 
                        
                            % ptojected distance 
                            dp = abs(ep(j) - pep(k));
                            % actual distance 
                            dist = (abs(dp*cosd(ang_j))+abs(dp*cosd(ang_k)))/2;
                            if(dist >= min_dist)
                                CL.ace{d}(round((ep(j)+pep(k))/2) ,i) = 1;
                                CL.AparEdgeLocX{d}(round((ep(j)+pep(k))/2) ,i, :) = [ep(j) i];
                                CL.AparEdgeLocY{d}(round((ep(j)+pep(k))/2) ,i, :) = [pep(k) i];
                            end
                    end

                end
            end
        end
    end
    
    I_ace = CL.ace{1} | CL.ace{2}';    
   
    if (cc_th ~= 0)
        Iout = cc_threshold(I_ace, cc_th, 1);
    else
        Iout = I_ace;
    end
    
    if(display == 1) 
        imwrite(1 - CL.ace{1},       sprintf('ACE/Results_2013_02_02/%s_2_ace_apar_edges_horizontal.png', img_name));
        imwrite(1 - CL.ace{2}',      sprintf('ACE/Results_2013_02_02/%s_3_ace_apar_edges_vertcal.png', img_name));
        imwrite(1 - I_ace,           sprintf('ACE/Results_2013_02_02/%s_4_ace_apar_edges_combined.png', img_name));
        imwrite(1 - Iout,            sprintf('ACE/Results_2013_02_02/%s_5_ace_apar_edges_refined.png', img_name));
        imwrite(1 - CL.CannyEdge{1}, sprintf('ACE/Results_2013_02_02/%s_1_ace_edge_map.png', img_name));
    end
    
%     
%     image = zeros(imH, imW, 3);
%     image(:,:,2) = img;
%     image(:,:,1) = I_ace; 
%     
%     figure;imshow(image,[]);
%     
end


function ang = GetAngle(x, y)

    if x<0 && y<0
        ang = 180 + atand(double(abs(y)/abs(x)));
    elseif x<0 && y>=0
        ang = 180 - atand(double(    y /abs(x)));
    elseif x>=0 && y<0
        ang = 360 - atand(double(abs(y)/    x));
    else
        ang =       atand(double(    y /    x));
    end    
end


function Ang = GetAngleM(X, Y)
        
   Ang = 180 * ((X<0)&(Y< 0) | ((X< 0)&(Y>=0))) + 360*((X>=0)&(Y<0)) ...
             + ((X<0)&(Y< 0) | ((X>=0)&(Y>=0))) .* atand(double(abs(Y)./abs(X))) ...
             - ((X<0)&(Y>=0) | ((X>=0)&(Y< 0))) .* atand(double(abs(Y)./abs(X)));

end

