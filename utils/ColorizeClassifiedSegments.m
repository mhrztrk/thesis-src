function [col_image] = ColorizeClassifiedSegments(original_image, segment_labels, label_classes)
    
    col_image = original_image;
    rgn_cnt = max(max(segment_labels));
    imH = size(original_image,1);
    imW = size(original_image,2);

    for k = 1:rgn_cnt
        % check whether the segment is labelled as road
        if (label_classes(k) == 1) 
            % If so, generate segment mask
            msk = (segment_labels == k);
            % According to this mask, colorize the segment
            % Colorization means removing Green and Blue components of that
            % segment, this make segment looks red.
            
            for i = 1:imH
                for j = 1:imW
                    if (msk(i,j) == 1)
                        col_image(i,j,1) = original_image(i,j,1);
                        col_image(i,j,2) = 0;
                        col_image(i,j,3) = 0;
                    end
                end
            end
            
%             [indx,indy] = find(segment_labels==k);
%             col_image(indx,indy,1) = original_image(indx,indy,1);
%             col_image(indx,indy,2) = 0;
%             col_image(indx,indy,3) = 0;
            
        end
    end

end