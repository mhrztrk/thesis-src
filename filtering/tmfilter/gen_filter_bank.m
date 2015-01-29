%
% ws - window size 
% min_rw - minimum road width
% max_rw - maximum road width
%
% both must be odd number & ws > 3*max_rw
%
function fb = gen_filter_bank(ws, min_rw, max_rw)

    fb = zeros(ws, ws, ((max_rw-min_rw)/2+1));
    gf = fspecial('gaussian', [ws ws], ws/3);

    for i=min_rw:2:max_rw
        x = zeros(1,ws);
        x(((ws-3*i)/2+1):((ws+3*i)/2)) = 1:(3*i);
        y = (-1 * sin(pi*x/i))*(max_rw/i);
        %y = (-1 * sin(pi*x/i));
        fb(:,:,(i-min_rw)/2+1) = repmat(y,ws,1).*gf;
    end

end
