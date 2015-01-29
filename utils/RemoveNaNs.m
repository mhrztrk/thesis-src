function [out] = RemoveNaNs(in,dim)
    msk = isnan(in(:,dim));
    m = find(msk == 0);
    out = in(m,:);
end