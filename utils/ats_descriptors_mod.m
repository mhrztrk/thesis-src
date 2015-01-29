% 
% Automated Road Network Extraction from High Spatial Resolution Multi-Spectral Imagery (Thesis)
% URL:http://www.geomatics.ucalgary.ca/research/publications/GradTheses.html
%
function [mean comp ecc var] = ats_descriptors_mod(I, w, step)

    % get image size
    [imH imW] = size(I);

    % Initialize variables
    mean = zeros(imH, imW);
    comp = zeros(imH, imW);
    ecc  = zeros(imH, imW);
    var  = zeros(imH, imW);
    
    % w must be odd, if not, just make it odd
    if(mod(w,2)~=1)
        w = w + 1;
    end

    % Create road filter    
    f = gen_filter_bank(3*w,w,w);
    fhsz = (size(f,1)-1)/2;
    
    % Compute all rotated masks
    ang = 360 / step;
    rmasks = zeros(size(f,1), size(f,2), step);
    for i=1:step
        rmasks(:,:,i) = imrotate(f, ang*i, 'crop');
    end
    
    for i=(fhsz+1):(imH-fhsz)
        for j=(fhsz+1):(imW-fhsz)
            pval = zeros(1,step);
            for k=1:step
                B = I((i-fhsz):(i+fhsz),(j-fhsz):(j+fhsz)).*rmasks(:,:,k);
                pval(k) = sum(B(:));
            end
            [X Y] = GetPolygon(pval);
            mean(i,j) = CalcMean(pval);
            comp(i,j) = CalcComp(X, Y);
            ecc(i,j)  = CalcEcc(X, Y);
            var(i,j)  = CalcVar(pval);
        end
    end
   
end


function [X Y] = GetPolygon(data)
    ang = 360 / size(data,2);
    X = zeros(1, size(data,2));
    Y = X;
    stp = size(data,2);
    
    for i=1:stp
        X(i) = data(i) * cos((ang * i)*pi/180);
        Y(i) = data(i) * sin((ang * i)*pi/180);
    end
end


function [M] = CalcMean(data)
    M = mean(data);
end

function [V] = CalcVar(data)
    M = CalcMean(data);
    V = sum((data - M).*(data - M));
end


function [C] = CalcComp(X,Y)

    stp = size(X,2);

    % calculate ares of polygon
    A=polyarea(X,Y);
    
    % calculate perimeter of polygon
    P = 0;
    for i=1:(stp-1)
        P = P + sqrt((X(i)-X(i+1))^2+(Y(i)-Y(i+1))^2);
    end
    P = P + sqrt((X(stp)-X(1))^2+(Y(stp)-Y(1))^2);
    
    % calculate compactness
    C = (4*pi*A)/(P^2);
    
%     plot(X,Y);grid;hold on;
%     scatter(0,0,'r','x');
%     scatter(X,Y,'r','x');hold off;
%     drawnow;
%     waitforbuttonpress;
end

function [E] = CalcEcc(X, Y)
    mX = mean(X);
    mY = mean(Y);
    E = sqrt(mX^2+mY^2);
end