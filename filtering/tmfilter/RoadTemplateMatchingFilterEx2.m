%   
%   I :  binary or gray-level image
%   N :  iteration count
%
%	road_range : [min_road_width max_road_width] (value of min & max must be be odd numbers )
% 	show_results : 	0 -> display nothing
%					1 -> display result of each iteration
%					2 -> save each result of each iteration 
%
function [result] = RoadTemplateMatchingFilterEx2(I, N, angle_step, road_range, show_results, method)
   
	% image name 
    filename = 'img1_tmf';

	% generate filter bank
    ws = 3 * road_range(2);
    fb = gen_filter_bank_int2(ws, road_range(1), road_range(2));

	% normalize image
    Isc = mat2gray(I);
    
    rot_step_cnt = 180 / angle_step + 1;

	% window half size
    wsn = (ws+1)/2;
    
	%% generate kernels
    kernel = zeros(ws,ws,size(fb,3),rot_step_cnt);   
    kernel_left = zeros(ws,ws,size(fb,3),rot_step_cnt);   
    kernel_right = zeros(ws,ws,size(fb,3),rot_step_cnt);   
    
    for m=1:size(fb,3)
        for i=1:rot_step_cnt
            kernel(:,:,m,i) = imrotate(fb(:,:,m),(i-1)*angle_step,'crop');
            
            kernel_left(:,1:wsn,m,i) = fb(:,1:wsn,m);
            kernel_left(:,:,m,i) = imrotate(kernel_left(:,:,m,i),(i-1)*angle_step,'crop');
            
            kernel_right(:,wsn:end,m,i) = fb(:,wsn:end,m);
            kernel_right(:,:,m,i) = imrotate(kernel_right(:,:,m,i),(i-1)*angle_step,'crop');

        end
    end
	%%
	
    if (method == 2) || (method == 3)
        kernel_right = filter_balance(kernel_right);
        kernel_left = filter_balance(kernel_left);
    end
    
    result = Isc;
	
    for k=1:N         
        imfilt = -Inf*ones(size(result,1),size(result,2),2);
        
        kk = max(4 - k,1);
        
        %for m=kk:(size(fb,3) + kk - 3)  
        for m=1:(size(fb,3) + kk - 3) 
    
            for i=1:rot_step_cnt
                 
                if (method == 1)
                    
                    imfilt(:,:,1) = imfilter(result, kernel(:,:,m,i));
                    
                elseif (method == 2)
                    
                    % left response
                    res1 = imfilter(result, kernel_left(:,:,m,i));

                    % right response
                    res2 = imfilter(result, kernel_right(:,:,m,i));     
                    
                    imfilt(:,:,1) = min(res1,res2);
                    
                elseif (method == 3)
                    
                    % left response
                    res1 = imfilter(result, kernel_left(:,:,m,i));

                    % right response
                    res2 = imfilter(result, kernel_right(:,:,m,i));
                
                    res = imfilter(result, kernel(:,:,m,i));
                    
                    imfilt(:,:,1) = max(min(res1,res2),res);
                
                end

                in_result = max(imfilt,[],3);
                imfilt(:,:,2) = in_result;
            end
        end
        
        result = mat2gray(in_result);

        if(show_results == 1)
            figure;imshow(result,[]);
        elseif (show_results == 2)
            imwrite(mat2gray(result), sprintf('%s_iter%d.png', filename,k));
        end
        
    end
    
end


%
% ws - window size 
% min_rw - minimum road width
% max_rw - maximum road width
%
% both min_rw and max_rw must be odd number & ws >= 3*max_rw
%
function fb = gen_filter_bank_int(ws, min_rw, max_rw)

    fb = zeros(ws, ws, ((max_rw-min_rw)/2+1));
    gf = fspecial('gaussian', [ws ws], ws/3);
    
    for i=min_rw:2:max_rw
        x = zeros(1,ws);
        x(((ws-3*(i-1))/2+1):((ws+3*(i-1))/2)) = 1:(3*(i-1));
        %y = (-1 * sin(pi*x/i));
        y = (-1 * sin(pi*x/(i-1)))*(max_rw/i);
        %y = (-1 * sin(pi*x/(i-1))); y = y ./ sum(y(y>0));
        fb(:,:,(i-min_rw)/2+1) = repmat(y,ws,1).*gf;
    end

end



function fb = gen_filter_bank_int2(ws, min_rw, max_rw)

    fb = zeros(ws, ws, ((max_rw-min_rw)/2+1));
    gf = fspecial('gaussian', [ws ws], ws/3);
    
    x_n = 0:4;
    y_n = -1 * sin(pi*x_n/4);
    y_n = y_n ./ (2*sum(abs(y_n)));
    
    for i=min_rw:2:max_rw
        x = zeros(1,ws);
        ws2 = (ws + 1)/2;
        
        x((ws2-(i-1)/2):(ws2+(i-1)/2)) = (i-1):(2*(i-1));
        
        y = (-1 * sin(pi*x/(i-1))); y = y ./ sum(y(y>0));
        y((ws2-(i+7)/2):(ws2-(i-1)/2)) = y_n;
        y((ws2+(i-1)/2):(ws2+(i+7)/2)) = y_n;
        
        fb(:,:,(i-min_rw)/2+1) = repmat(y,ws,1).*gf;
    end

end



function [filtb] = filter_balance(filt)

    fsz = size(filt);
    filtb = zeros(fsz);
    
    for i=1:fsz(3)
        for j=1:fsz(4)
            krn = filt(:,:,i,j);
            filtb(:,:,i,j) = (1/sum(krn(krn>0)))*(krn.*(krn>0)) + ...
                          abs(1/sum(krn(krn<0)))*(krn.*(krn<0));
                        
        end
    end

end

