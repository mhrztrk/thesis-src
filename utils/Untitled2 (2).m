 ks = size(filtb);
 
 Sp = zeros(ks(3),ks(4));
 Sn = zeros(ks(3),ks(4));
 
for i=1:ks(3)
    for j=1:ks(4)
        krn = filtb(:,:,i,j);
        
        krn_p = krn(krn>0);
        krn_n = krn(krn<0);
       
        Sp(i,j) = sum(krn_p);
        Sn(i,j) = abs(sum(krn_n));
        
    end
end