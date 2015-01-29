function [] = MRF()
    % im = ???
    min_dist = 8;   % minimum distance between parallel edges.
    max_dist = 25;  % maximim distance between parallel edges
    min_defl = 20;  % max angle deflection between apar edge points' gradient directions 

    % Apply ACE algorithm to find (possible) road centerlines 
    aparedges = ace( im, min_dist, max_dist, min_defl);

    % Link edge pixels together into lists of sequential edge points, one
    % list for each edge contour. A contour/edgelist starts/stops at an 
    % ending or a junction with another contour/edgelist.
    % Here we discard contours less than 10 pixels long.
    [edgelist, labelededgeim] = edgelink(aparedges, 10);

    % Fit line segments to the edgelists
    tol = 2;  % Line segments are fitted with maximum deviation from
              % original edge of 2 pixels.
    seglist = lineseg(edgelist, tol);

    k = 1;
    for i=1:size(seglist,2)
        for j=1:(size(seglist{i},1)-1)
            X(k) = (seglist{i}(j,1) + seglist{i}(j+1,1))/2;
            Y(k) = (seglist{i}(j,2) + seglist{i}(j+1,2))/2;
            k = k + 1;
        end
    end


    % Create Adjacency Graph


    % find maximal cliques
    Cliques = maximalCliques(Adj);

    nconn = size(Adj,1);    % # of connections
    
    deltaE = Inf;
    
    step = 1;
    c = 0.1;
    T = 100000;
    t = -Inf;
    K = 0;
    
    labels = randi([0 1], [nconn 1]);
    
    while(deltaE > t && K <= 20)
        
        for k=1:step
            for i=1:nconn

                labels_1 = labels;
                labels_2 = labels;

                labels_1(i) = 0;
                labels_2(i) = 1;

                E_1 = exp(-CalculateEnergy(Segments, labels_1, Cliques, D, mrf)/T);
                E_2 = exp(-CalculateEnergy(Segments, labels_2, Cliques, D, mrf)/T);

                p = E_1 / (E_1 + E_2);

                r = randi([0 1]);

                if(r < p)
                    labels = labels_1;
                    E_new = E_1;
                else
                    labels = labels_2;
                    E_new = E_2;
                end

                E_old = E_new;
                EE(k) = E_new;
                
                % Calculate total energy
                % E_new = exp((sum(V) + sum(Vc))/T);

                deltaE = abs(E_new - E_old);    % stop when energy change is small
            end
        end
        T = T * c;      % Decrease temperature
        K = K + 1;      % Advance iteration counter
        
    end
    
end
