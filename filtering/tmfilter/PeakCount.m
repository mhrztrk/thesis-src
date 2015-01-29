function [pcmap] = PeakCount(I, Template, Th)

    [H W] = size(I);
    tw = size(Template,1); 

    angle_step   = 15;
    rot_step_cnt = 360 / angle_step; 
    
    pcmap = zeros(H,W);
    
    kernel = zeros(tw,tw,rot_step_cnt);
     for i=1:rot_step_cnt
        kernel(:,:,i) = imrotate(Template,i*angle_step,'crop');
     end
    
    for m=((tw+1)/2):(H-((tw-1)/2))
        for n=((tw+1)/2):(W-((tw-1)/2))
            
            value = zeros(1, rot_step_cnt);
            for i=1:rot_step_cnt
                value(1, i) = ApplyFilter(I, kernel(:,:,i), n, m);
            end

            [vals locs]  = FindMaximaLocs(value, 1, 180, rot_step_cnt);
            pcmap(m,n) = size(locs,2);
            
        end
    end
    


end


%%
function [fout] = ApplyFilter(img, filter, pt_x, pt_y)

        nim = uint16((size(filter,1)-1)/2);
        
%       fout = conv2(img((pt_y - nim):(pt_y + nim),...
%                          (pt_x - nim):(pt_x + nim)), ...
%                             filter,'valid');
%       value(1,1:end) = value(1,linspace(rot_step_cnt,1,rot_step_cnt));

        fout = sum(sum(img( (pt_y - nim):(pt_y + nim),...
                            (pt_x - nim):(pt_x + nim)).*filter));
                                         
end