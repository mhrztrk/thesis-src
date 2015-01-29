function [Rmin, qmetrics] = ocsvm_classification(img, mask, ref, imgName, baseDir)

	groups = ones(size(Y,1), 1);
	nSample = size(groups,1);
	
	
	% cross-validation
	
	bestnu = param(1);  %% nu is critical, change if feature set is changed.
	bestg = 16;
	
	bestcv = 0;
	i = 1;
	
	for log2g = -1:10
		
		qmetrics{i}.g = 2^log2g;
		qmetrics{i}.nu = bestnu;
		
		
		cmd = ['-q -v 5 -s 2 -n ', num2str(bestnu), ' -g ', num2str(2^log2g)];
		cv = svmtrain(double(groups), double(Y), cmd);
		if (cv >= bestcv),
			bestcv = cv; bestg = 2^log2g;
		end

		cmd = ['-s 2 -t 2 -n ', num2str(bestnu), ' -g ', num2str(2^log2g)];
		models = svmtrain(double(groups), double(Y), cmd);

		fprintf('total SVs = %d\n', models.totalSV);
		fprintf('%g %g (best nu=%g, g=%g, rate=%g)\n', log2g, cv, bestnu, bestg, bestcv);
	
		[predicted_label, ~, ~] = svmpredict(double(ones(size(Y,1),1)), double(Y), models);
		
		% number of outliers
		nOut = sum(predicted_label <= 0);
		
		[~, acc, decision_values] = svmpredict(double(ones(size(X,1),1)), double(X), models);
		
		R = reshape(decision_values, size(I,1), size(I,2));
		
		[qmetrics{i}.completeness, qmetrics{i}.correctness] = fPreRecallHesapla_v3(ref, R>0 , '', 0);
		
		qmetrics{i}.quality = (qmetrics{i}.completeness*qmetrics{i}.correctness)/...
			(qmetrics{i}.completeness - qmetrics{i}.completeness*qmetrics{i}.correctness + qmetrics{i}.correctness);
		
		strPerEval = sprintf('Completeness = %.4f, Correctness = %.4f, Quality = %.4f', ...
			qmetrics{i}.completeness, qmetrics{i}.correctness, qmetrics{i}.quality);
		
		strTitle = sprintf('@(nu=%5.4f, g=%6.1f) => cvAcc=%4.1f , clAcc=%4.1f , mse=%4.2f, fsv = (%g/%g), outliers=(%d/%d)',...
			bestnu, 2^log2g, cv, acc(1), acc(2), models.totalSV, size(Y,1),nOut,nSample);
		
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

	cmd = ['-s 2 -t 2 -n ', num2str(bestnu), ' -g ', num2str(bestg)];
	models = svmtrain(double(groups), double(Y), cmd);

	[~, ~, decision_values] = svmpredict(double(ones(size(X,1),1)), double(X), models);
	
	R = reshape(decision_values, size(I,1), size(I,2));

	
	if(size(Y,2) == 2)
		distsurfplot(Y, models,range, R>0, '');
	end
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
