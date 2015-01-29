function [linelist newLineNr] = DivideLineIntoTwo(linelist, lineNr, cutPt)

    nline = size(linelist, 2);
    
    nline = nline + 1;
    
    linelist(nline).nr = nline;
    linelist(nline).s  = cutPt;
    linelist(nline).e  = linelist(lineNr).e;
    linelist(nline).c  = (linelist(nline).s + linelist(nline).e)/2;   
     
    sp_conns = linelist(lineNr).conn(linelist(lineNr).connpt == 0);
    ep_conns = linelist(lineNr).conn(linelist(lineNr).connpt == 1);

    for i=1:size(ep_conns,2)
        linelist(ep_conns(i)).conn(linelist(ep_conns(i)).conn == lineNr) = nline;
    end
    
    linelist(nline).conn   = ep_conns;
    linelist(nline).nconn  = size(ep_conns,2);
    linelist(nline).connpt = ones(1,linelist(nline).nconn);
    
    linelist(nline).nconn  = linelist(nline).nconn + 1;
    linelist(nline).conn(linelist(nline).nconn) = lineNr;
    linelist(nline).connpt(linelist(nline).nconn) = 0;
    
    linelist(lineNr).e = cutPt;
    linelist(lineNr).c = (linelist(lineNr).s + linelist(lineNr).e)/2;
    
    linelist(lineNr).conn   = sp_conns;
    linelist(lineNr).nconn  = size(sp_conns,2);
    linelist(lineNr).connpt = zeros(1,linelist(lineNr).nconn);
    
    linelist(lineNr).nconn  = linelist(lineNr).nconn + 1;
    linelist(lineNr).conn(linelist(lineNr).nconn) = nline;   
    linelist(lineNr).connpt(linelist(lineNr).nconn) = 1;
    
    newLineNr = lineNr;
    
end