function [result] = RoadTemplateMatchingFilter(bin_image, pass_count, template, angle_step, show_results)
   
    bin_image = mat2gray(bin_image);
    figure;imshow(bin_image,[]);
    
    rot_step_cnt = 180 / angle_step + 1;
    result = bin_image;
    
    szTemplate = size(template,1);
    kernel = zeros(szTemplate,szTemplate,rot_step_cnt);
    
    for i=1:rot_step_cnt
        kernel(:,:,i) = imrotate(template,(i-1)*angle_step,'crop');
    end
   
    for k=1:pass_count         
        imfilt = -Inf*ones(size(result,1),size(result,2),2);
        for i=1:rot_step_cnt
            imfilt(:,:,1) = imfilter(result, kernel(:,:,i));
            in_result = max(imfilt,[],3);
            imfilt(:,:,2) = in_result;
        end

        %result = max(imfilt,[],3);
        result = mat2gray(in_result);
        
        if(show_results == 1)
            figure;imshow(result,[]);
        elseif (show_results == 2)
            imwrite(mat2gray(result), sprintf('%s_iter%d.png', filename,k));
        end
        
    end
    
end

