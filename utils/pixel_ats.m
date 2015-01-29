function [A] = pixel_ats(I, p, f, as)

    rs_cnt = 180 / as; 
    A = zeros(1, rs_cnt);
     
     for i=1:rs_cnt
        krnl = imrotate(f, i*as,'crop');
        A(1, i) = ApplyFilter(I, krnl, p(1), p(2));
     end
     
end

function [fout] = ApplyFilter(img, filter, x, y)
        fhsz = uint16((size(filter,1)-1)/2);
        img_wnd = img((x-fhsz):(x+fhsz),(y-fhsz):(y+fhsz)); 
        fr = imfilter(img_wnd, filter);
        fout = fr(fhsz+1,fhsz+1);
end