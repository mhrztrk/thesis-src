%   
%   I :  binary or gray-level image
%   N :  iteration count
%
%
function [result] = RoadTemplateMatchingFilterEx(I, N, angle_step, road_range, show_results)
   
    ws = 3 * road_range(2);
    fb = gen_filter_bank_int2(ws, road_range(1), road_range(2));

    Isc = mat2gray(I);
    figure;imshow(Isc,[]);
    
    rot_step_cnt = 180 / angle_step + 1;

    kernel = zeros(ws,ws,size(fb,3),rot_step_cnt);   
    for m=1:size(fb,3)
        for i=1:rot_step_cnt
            kernel(:,:,m,i) = imrotate(fb(:,:,m),(i-1)*angle_step,'crop');
        end
    end

    result = Isc;
    for k=1:N         
        imfilt = -Inf*ones(size(result,1),size(result,2),2);
        for m=1:size(fb,3)
            for i=1:rot_step_cnt
                imfilt(:,:,1) = imfilter(result, kernel(:,:,m,i));
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



