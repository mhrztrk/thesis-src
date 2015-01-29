% 
% I : input image
% T : mask
% method : 
%
function [clres qmetrics] = imclassify(I, M, ref, method, param, range, imgName, baseDir)

    % image size
    [imH imW imD] = size(I);
    
    X = reshape(I, imH*imW, imD);   % image pixels in vector form
    
    % transformation 
    % X = transform_to_sth(X);
    
    Y = X((reshape(M, imH*imW, 1) == 1),:); % training set
    
%     fid = fopen('training_set_file','w+');
%     
%     for i=1:size(Y,1)
%         fprintf(fid, '+1 1:%f 2:%f\r\n', Y(i,1), Y(i,2));
%     end
    
    % remove NaN elements
    NanMask = false(size(Y,1),1);
    for i=1:imD
        NanMask = NanMask | isnan(Y(:,i));
    end
    
    for i=1:imD
        tmp = Y(:,i);
        Z(:,i)= tmp(~NanMask);
    end
    
    Y = Z;
    
    outDir = sprintf('%s/%s', baseDir, imgName);
    mkdir(outDir);
    fid = fopen(sprintf('%s/results_%s.txt', outDir, imgName), 'a');
    
    if strcmp(method,'mahal')
        d = mahal(double(X), double(Y)); % Mahalanobis
        D = 1 - (d - min(d))/(max(d) - min(d));
        R = reshape(D, imH, imW);
        
    elseif strcmp(method,'gmm')
        %options = statset('Display','final');
        %obj = gmdistribution.fit(Y,2,'Options',options);
         
        obj = gmdistribution.fit(Y, param);
        
        if(param ~= 1)
            [P] = posterior(obj,X);

            R = zeros(imH, imW, param);
            
            fprintf('Component ratios: ')
            
            for i=1:param
                R(:,:,i) = reshape(P(:,i), size(I,1), size(I,2));
                fprintf('%.3f ', obj.PComponents(i));
            end
            fprintf('\n');
            
        else
            R = pdf(obj,X);
            R = reshape(R,1000,1000);
        end
        
        if(size(Y,2) == 2)
            gmmcontour(Y, obj);
        end
        
        
        
    elseif strcmp(method,'bhattacharya')
        R = BhattacharyyaClassifier(I, M, 5);
		
    elseif strcmp(method,'ocsvm')

        groups = ones(size(Y,1), 1);
        nSample = size(groups,1);
                
        % cross-validation
        
        bestnu = param(1);  %% nu is critical, change if feature set is changed.
        
        bestcv = 0;
        i = 1;
        
        cv = zeros(1,10);
        clres = [];
        
        for log2g = -1:10
            
            qmetrics.g = 2^log2g;
            qmetrics.nu = bestnu;

            cmd = ['-q -v 5 -s 2 -n ', num2str(bestnu), ' -g ', num2str(2^log2g)];
            cv(i) = svmtrain(double(groups), double(Y), cmd);

            % if reduction in cv is more than 10%, stop searching 
            if((cv(1) - cv(i)) > 0.05*cv(1))
                clres = (R > 0);
                qmetrics.g = 2^(log2g-1);
                break;
            end
            
            cmd = ['-s 2 -t 2 -n ', num2str(bestnu), ' -g ', num2str(2^log2g)];
            models = svmtrain(double(groups), double(Y), cmd);

            fprintf('total SVs = %d\n', models.totalSV);
        
            [predicted_label, ~, ~] = svmpredict(double(ones(size(Y,1),1)), double(Y), models);
            
            % number of outliers
            nOut = sum(predicted_label <= 0);
            
            [~, acc, decision_values] = svmpredict(double(ones(size(X,1),1)), double(X), models);
            
            R = reshape(decision_values, size(I,1), size(I,2));
            
            [qmetrics.completeness, qmetrics.correctness] = fPreRecallHesapla_v3(ref, R>0 , '', 0);
            
            qmetrics.quality = (qmetrics.completeness*qmetrics.correctness)/...
                (qmetrics.completeness - qmetrics.completeness*qmetrics.correctness + qmetrics.correctness);
            
            strPerEval = sprintf('Completeness = %.4f, Correctness = %.4f, Quality = %.4f', ...
                qmetrics.completeness, qmetrics.correctness, qmetrics.quality);
            
            strTitle = sprintf('@(nu=%5.4f, g=%6.1f) => cvAcc=%4.1f , clAcc=%4.1f , mse=%4.2f, fsv = (%g/%g), outliers=(%d/%d)',...
                bestnu, 2^log2g, cv(i), acc(1), acc(2), models.totalSV, size(Y,1),nOut,nSample);
            
            if(size(Y,2) == 2)
                [img fds] = distsurfplot(Y, models, range, R>0, strTitle);
                imwrite(img, sprintf('%s/%2d - ocsvm_decboundary_%s_(nu_%g)_(g_%g).png', outDir, i, imgName, bestnu, 2^log2g));
            else
                fds = 0;
            end
            
            fprintf(fid, '%2d - %s, fsd = %4.1f, %s\n', i, strTitle, fds, strPerEval);
            
            imwrite(R>0, sprintf('%s/%2d - ocsvm_clres_%s_(nu_%g)_(g_%g).png', outDir, i, imgName, bestnu, 2^log2g));
            
            i = i + 1;
            
        end

        if(isempty(clres))
            clres = R > 0;
        end
        
        figure;plot(cv);
        
        cmd = ['-s 2 -t 2 -n ', num2str(qmetrics.nu), ' -g ', num2str(qmetrics.g)];
        models = svmtrain(double(groups), double(Y), cmd);

        [~, ~, decision_values] = svmpredict(double(ones(size(X,1),1)), double(X), models);
        
        R = reshape(decision_values, size(I,1), size(I,2));

        if(size(Y,2) == 2)
            distsurfplot(Y, models,range, R>0, '');
        end
        
    end

    fclose(fid);
    
end




function gmmcontour(data, gmm)

    x1 = 0.001:0.001:1;
    x2 = 0.001:0.001:1;

    [X1 X2] = meshgrid(x1,x2);

    figure;scatter(data(:,1), data(:,2), 10, 'o', 'b', ...
            'MarkerFaceColor',[0.2 0 1],...
            'MarkerEdgeColor',[0.2 0 1]); 
    
    hold on
        
    for i=1:gmm.NComponents;
        O{i} = pdf(gmdistribution(gmm.mu(i,:), gmm.Sigma(:,:,i)),[X1(:) X2(:)]);
        contour(x1,x1,reshape(O{i},size(x1,2),size(x1,2)));
    end
           
    % Create xlabel
    xlabel({'','Principal Component - 1'});

    % Create ylabel
    ylabel({'Principal Component - 2',''});

end



function [img, FDS] = distsurfplot(data, models,range, result, strTitle)

    %%

     x1 = linspace(range(1),range(2),1000);
     x2 = x1;
     
     d = (range(2)-range(1))/1000;
     
%%
    pcnt = size(x1,2);
     
    [X1 X2] = meshgrid(x1,x2);

    [~, ~, decision_values] = svmpredict(double(ones(pcnt*pcnt,1)), [X1(:) X2(:)], models);
 
    F = reshape(decision_values,pcnt,pcnt);
  
    % fraction of decision space
    FDS = round(sum(decision_values > 0) * 1000 / size(decision_values,1))/10;
    
    f = figure;
    
    subplot(1,2,2);imshow(result, []);
    title(strTitle);
    
    screen_size = get(0, 'ScreenSize');

    set(f, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    
    img = 0.5 * (1 - (F>0)) + (1-bwmorph(edge(F>0,'canny'),'thicken'));
    
    subplot(1,2,1);imshow(img,[]);
    
    hold on; scatter(   (data(:,1) - range(1))/d, ...
                        (data(:,2) - range(1))/d, ...
                        10, 'x', 'b','MarkerFaceColor',[0.2 0 1],...
                        'MarkerEdgeColor',[0.2 0 1]);  
    
    % contour(x1,x2,F, 'ShowText','on','LineWidth',2,'Fill','on');
    
    %colormap('default');
        
    
    f = figure('units','normalized','outerposition',[0 0 1 1]);
        
    imshow(img,[]);
      
    hold on; scatter(   (data(:,1) - range(1))/d, ...
                        (data(:,2) - range(1))/d, ...
                        10, 'x', 'b','MarkerFaceColor',[0.2 0 1],...
                        'MarkerEdgeColor',[0.2 0 1]);       
    
    grid;
                    
    title(strTitle);
    
    set(f, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    
    tmp = getframe(gcf);
    img = tmp.cdata;
    
    close 

    %%
end


function [A] = decisionarea(obj, dist, d)
    
    imH = d(1);
    imW = d(2);
    
    x1 = 0.001:0.001:1;
    x2 = 0.001:0.001:1;

    [X1 X2] = meshgrid(x1,x2);

    [P] = posterior(obj, [X1(:) X2(:)]);

    for i=1:obj.NComponents
        R{i} = reshape(P(:,i), imH, imW);
    end    
    
    for i=1:obj.NComponents
        O{i} = reshape(pdf(gmdistribution(obj.mu(i,:), obj.Sigma(:,:,i)),[X1(:) X2(:)]),...
            size(x1,2),size(x1,2));
    end
    
    dec = ones(size(x1,2),size(x1,2),obj.NComponents);
    
    for i=1:obj.NComponents
        j=1:obj.NComponents;
        j = j(j ~= i);
        for k=j
            dec(:,:,i) = dec(:,:,i) & (O{i} > O{k}); 
        end  
    end
    
    decb = zeros(size(x1,2),size(x1,2));
    for i=1:size(dist.comp,2)
        decb = decb | dec(:,:,dist.comp(i));
    end
    
    % figure;imshow(decb,[]);
    
    rp = regionprops(decb, 'area');
    A = sum([rp.Area]);  
    
end
