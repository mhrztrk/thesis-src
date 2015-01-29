%
% Pixel Shape Index
%
function [psi] = impsi(im,th)
    
    % Spectral homogenity threshold.
    T1 = th; 

    % Total # of pixels threshold in the direction line.
    T2 = 15;
    
    % Image size
    [imH imW imS] = size(im);
 
    % Angle step between direction lines.
    ang_step = 30;
    step_cnt = 360 / ang_step;
    
    psi = zeros(imH,imW);
    
    % pixel locations in the line 
    lines = zeros(step_cnt, T2, 2);
    
    lin = 1:T2;
    
    for k=1:step_cnt
        lines(k,:,1) = round(lin * sin((ang_step*(k-1))*pi/180));
        lines(k,:,2) = round(lin * cos((ang_step*(k-1))*pi/180));
    end
    
    L = zeros(step_cnt,1);
    
    for i=(T2+1):(imH-T2)
        for j=(T2+1):(imW-T2)
            for k=1:step_cnt
                for m=1:T2
                    ph = sum((im(lines(k,m,1)+i,lines(k,m,2)+j,:)-im(i,j,:)).^2);
                    if(ph >= T1)
                        break;
                    end
                end

                % calculate length of this direction line.
                L(k) = sqrt(lines(k,m,1)^2 + lines(k,m,2)^2);
            end
            psi(i,j) = sum(L(:));
        end
    end
end