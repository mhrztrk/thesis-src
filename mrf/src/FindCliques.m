function [Cliques] = FindCliques(linelist)
   
    ncliq = 0;
   
    nline = size(linelist,2);
    
    vmat = zeros(nline, 2); 
    
    
    for i=1:nline
         
        if(ncliq == 53)
            disp('a');
        end
        
        if(vmat(i,1) == 0)
            if(linelist(i).nconn)
                sp_conns = linelist(i).conn(linelist(i).connpt == 0); 
                if(~isempty(sp_conns))
                    for k=1:size(sp_conns,2)
                        ind = find(linelist(sp_conns(k)).conn == i);
                        if(linelist(sp_conns(k)).connpt(ind) == 0)
                            vmat(sp_conns(k),1) = 1;
                        else
                            vmat(sp_conns(k),2) = 1;
                        end
                    end
                    ncliq = ncliq + 1;
                    Cliques{ncliq} = [i sp_conns]; 
                end
            else
                ncliq = ncliq + 1;
                Cliques{ncliq} = i;
            end
        end

        if(vmat(i,2) == 0)
            if(linelist(i).nconn)
                ep_conns = linelist(i).conn(linelist(i).connpt == 1);  
                if(~isempty(ep_conns))
                    for k=1:size(ep_conns,2)
                        ind = find(linelist(ep_conns(k)).conn == i);
                        if(linelist(ep_conns(k)).connpt(ind) == 0)
                            vmat(ep_conns(k),1) = 1;
                        else
                            vmat(ep_conns(k),2) = 1;
                        end
                    end
                    ncliq = ncliq + 1;
                    Cliques{ncliq} = [i ep_conns];                 
                end
            end
        end 
    end

    Cliques = Cliques';
end