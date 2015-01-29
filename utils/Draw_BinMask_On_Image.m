
function [MI] = Draw_BinMask_On_Image(I, M)

    MI = I(:,:,1:3);

    [X Y] = find(M == 1);
    MI(sub2ind(size(I), X, Y, 1*ones(size(X)))) = 1;
    MI(sub2ind(size(I), X, Y, 2*ones(size(X)))) = 0;
    MI(sub2ind(size(I), X, Y, 3*ones(size(X)))) = 0;

end