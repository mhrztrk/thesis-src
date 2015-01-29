function [result] = RoadTemplateMatchingFilter(bin_image, pass_count, template)

    %%
    %template = temp;
    %figure;imshow(template,[]);
    
    %pass_count = 5;
    %bin_image = (classified_image_2_1 + 1)/2;
    %filename = 'classified_image_2_1';
    
    figure;imshow(bin_image,[]);
    
    angle_step = 10;
    rot_step_cnt = 180 / angle_step + 1;
    result = bin_image;
    
    szTemplate = size(template,1);
    kernel = zeros(szTemplate,szTemplate,rot_step_cnt);
    
    for i=1:rot_step_cnt
        kernel(:,:,i) = imrotate(template,(i-1)*angle_step,'crop');
    end
    
    %%
    imfilt = zeros(size(result,1),size(result,2),rot_step_cnt);
    
    for k=1:pass_count         
        for i=1:rot_step_cnt
            imfilt(:,:,i) = imfilter(result, kernel(:,:,i));
        end
        
        result = max(imfilt,[],3);

        figure;imshow(result,[]);
        %imwrite(mat2gray(result), sprintf('%s_iter%d.png', filename,k));
    end

    %%
    clear imfilt;
    
    % remove the regions with negative value
%     for i=1:size(result,1)
%         for j = 1:size(result,2)
%             if result(i,j) < 0
%                 result(i,j) = 0;
%             end
%         end
%     end
%    
%     figure;imshow(result,[]);
end