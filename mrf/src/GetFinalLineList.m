function [fLineList] = GetFinalLineList(iLineList, labels)

    fLineList = iLineList(labels == 1);

    for i=1:size(fLineList,2)
        fLineList(i).connpt = fLineList(i).connpt(labels(fLineList(i).conn) == 1);
        fLineList(i).defl   = fLineList(i).defl(labels(fLineList(i).conn) == 1);
        fLineList(i).conn   = fLineList(i).conn(labels(fLineList(i).conn) == 1);
        fLineList(i).nconn  = size(fLineList(i).conn,2);
        fLineList(i).ncut = 0;
        fLineList(i).nclq = 0;
    end


    NrT = zeros(size(fLineList,2),1);

    for i=1:size(fLineList,2)
        NrT(i) = fLineList(i).nr;
        fLineList(i).nr = i;
    end
   
    for i=1:size(fLineList,2)
        for j=1:size(fLineList(i).conn,2)
            ind = find(NrT == fLineList(i).conn(j));
            fLineList(i).conn(j) = ind;
        end
        
    end
    
end