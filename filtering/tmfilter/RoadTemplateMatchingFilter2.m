function [result] = RoadTemplateMatchingFilter2(bin_image, pass_count, filter)

    %%
    template = filter;
    figure;imshow(template,[]);
    
    pass_count = 1;
    %bin_image = (classified_image_2_1 + 1)/2;
    %filename = 'classified_image_2_1';
    
    figure;imshow(bin_image,[]);
    
    angle_step = 15;
    rot_step_cnt = 360 / angle_step + 1;
    
    bend_angle_step = 15;
    min_bend_anle = 90;
    bending_step_cnt = (180 - min_bend_anle)/ bend_angle_step + 1;
    
    result = bin_image;
    
    %%
    for k=1:pass_count      
        
        imfilt2 = zeros(size(result,1),size(result,2),bending_step_cnt);
        
        for j=1:bending_step_cnt    
            
            bend_angle = min_bend_anle + (j-1)*bend_angle_step;
            
            tempA = imrotate(template, 90 - (bend_angle/2), 'crop');
            tempB = imrotate(template, - 90 + (bend_angle/2), 'crop');
            
            temp = template;
            temp(1:16,:) = tempA(1:16,:);
            temp(18:35,:) = tempB(18:35,:); 
            
            figure;imshow(temp,[]);
            
            imfilt = zeros(size(result,1),size(result,2),rot_step_cnt);

            for i=1:rot_step_cnt
                %template_rot = imrotate(template,(i-1)*angle_step,'crop');
                %kernel = template_rot(19:53,19:53);
                kernel = imrotate(temp,(i-1)*angle_step,'crop');
                imfilt(:,:,i) = imfilter(result, kernel);
            end
            
            imfilt2(:,:,j) = max(imfilt,[],3);
            
            clear imfilt;

            figure;imshow(imfilt2(:,:,j),[]);
            %imwrite(mat2gray(result), sprintf('%s_iter%d.png', filename,k));
        end
        
        result = max(imfilt2,[],3);
        figure;imshow(result,[]);         
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