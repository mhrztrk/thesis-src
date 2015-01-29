function [Rmin, qmetrics] = gmm_classification(img, mask, ref, imgName, baseDir)

	%% image size
	[imH imW imD] = size(img);

	X = reshape(img, imH*imW, imD);   % image pixels in vector form

	Y = X((reshape(mask, imH*imW, 1) == 1),:); % training set

	%%

	% minimum capture percentage, i.e at least climit percent must be captured by 
	% the mixture model.
	climit = 0.85;

	dist = [];

	outDir = sprintf('%s/%s', baseDir, imgName);
	mkdir(outDir);
	fid = fopen(sprintf('%s/results_%s.txt', outDir, imgName), 'a');

	nr_sol = 0;     % nr of solutions

	for kk=2:5
		
		obj = gmdistribution.fit(Y, kk);
	   
		fprintf('Component ratios: ')

		for i=1:kk
			fprintf('%.3f ', obj.PComponents(i));
		end
		fprintf('\n');    
		%%
		 
		n = 1;
		
		for i=1:(obj.NComponents-1)
		
			comb = combnk(1:obj.NComponents, i);
			
			for j=1:size(comb,1)
				totp = 0;
				for m=1:size(comb,2)
					totp = totp + obj.PComponents(comb(j,m));
				end
			
				if (totp > climit)

					tmp.success = 1;
					tmp.covered = round(totp*1000)/10;
					
					for m=1:size(comb,2)
						tmp.gauss{m} = gmdistribution(obj.mu(comb(j,m),:), obj.Sigma(:,:,comb(j,m)));
					end
					
					tmp.comp = comb(j,:);
					tmp.sparse = decisionarea(obj, tmp, [imH imW]);
					
					tmp.cva = CrossValidation(Y, obj, tmp.comp);
					
					if(tmp.sparse < 50.0)
						if(n==1)
							dist = tmp;
						else
							dist(n) = tmp;
						end
						n = n + 1; 
					end
					
				end
			end
		end    
		
	%     if(imD==2)
	%         plotgmm(obj, Y);
	%     end
		
		if(~isempty(dist))

			[P] = posterior(obj,X);

			R = [];

			skewmin = Inf;
			
			for i=1:kk
				R{i} = reshape(P(:,i), imH, imW);
			end    

			for i=1:size(dist,2)
				
				S = zeros(imH,imW);
				strComp = '[';
				for j=1:size(dist(i).comp, 2)
					S = S + R{dist(i).comp(j)};
					strComp = sprintf('%s %d', strComp, dist(i).comp(j));
				end
				strComp = sprintf('%s]', strComp);
				
				[Sdec] = classificationresult(dist(i), [imH imW], R);
				
				nr_sol = nr_sol + 1;
				
				[comp, corr] = fPreRecallHesapla_v3(ref, Sdec , '', 0);
				
				quality = comp*corr / (comp - comp*corr + corr);
				
				qmetrics{nr_sol}.ncomp      = [kk size(dist(i).comp, 2)];
				qmetrics{nr_sol}.decspace   = dist(i).sparse;
				qmetrics{nr_sol}.captured   = dist(i).covered;
				
				qmetrics{nr_sol}.completeness = comp;
				qmetrics{nr_sol}.correctness  = corr;
				qmetrics{nr_sol}.quality      = quality;
				
				
				figure; imshow(Sdec, []);
				
				fig_title = sprintf('N = %d, K = %d, cva = %4.1f, fds = %4.1f, captured = %4.1f, perf(Precision = %4.1f, Recall = %4.1f, Accuracy = %4.1f)', kk, size(dist(i).comp, 2),...
								dist(i).cva, dist(i).sparse, dist(i).covered,...
								Precision, Recall, Accuracy);
				title(fig_title);

				fprintf(fid, '%s\n', fig_title);
				
				if(imD == 2)
					[imDB] = drawdecisionboundary(obj, Y, dist(i), fig_title);
					 imwrite(imDB, sprintf('%s/%d_%d_decision_boundary.png', outDir, kk, i)); 
				end

				imwrite(Sdec, sprintf('%s/%d_%d_gmm_result.png', outDir, kk, i));
				
				if(sum(dist(i).sparse) < skewmin)
					Rmin = mat2gray(S);
					skewmin = sum(dist(i).sparse);
				end
			end
		end    
	end
	fclose(fid);
end


function A = ellipsearea(eig)
    A = pi*eig(1)*eig(2);
end


function plotgmm(obj, Y)

    x1 = 0.001:0.001:1;
    x2 = 0.001:0.001:1;

    [X1 X2] = meshgrid(x1,x2);
    
%     F1 = mvnpdf([X1(:) X2(:)],obj.mu(1,:),obj.Sigma(:,:,1));
%     F2 = mvnpdf([X1(:) X2(:)],obj.mu(2,:),obj.Sigma(:,:,2));
% 
%     F = reshape(F1,length(x2),length(x1)) + reshape(F2,length(x2),length(x1));
%     surf(x1,x2,F);

    for i=1:obj.NComponents
        O{i} = pdf(gmdistribution(obj.mu(i,:), obj.Sigma(:,:,i)),[X1(:) X2(:)]);
    end
    figure;scatter(Y(:,1), Y(:,2), 10, 'o', 'b',...
         'MarkerFaceColor',[0.2 0 1],...
        'MarkerEdgeColor',[0.2 0 1]); 

    hold on
    
    for i=1:obj.NComponents
        contour(x1,x1,reshape(O{i},size(x1,2),size(x1,2)));
    end
    
    % Create xlabel
    xlabel({'','Principal Component - 1'});

    % Create ylabel
    ylabel({'Principal Component - 2',''});

end


function [imDB] = drawdecisionboundary(obj, Y, dist, fig_title)

    x1 = 0.001:0.001:1;
    x2 = 0.001:0.001:1;

    
    [X1 X2] = meshgrid(x1,x2);

    for i=1:obj.NComponents
        O{i} = reshape(pdf(gmdistribution(obj.mu(i,:), obj.Sigma(:,:,i)),[X1(:) X2(:)]),...
            size(x1,2),size(x1,2));
    end
    
    figure;subplot(1,2,1);scatter(Y(:,1), Y(:,2), 10, 'o', 'b',...
         'MarkerFaceColor',[0.2 0 1],...
        'MarkerEdgeColor',[0.2 0 1]); 

    hold on
    
    for i=1:obj.NComponents
        contour(x1,x1,O{i});
    end

    dec = true(size(x1,2),size(x1,2),obj.NComponents);
    
    for i=1:obj.NComponents
        j=1:obj.NComponents;
        j = j(j ~= i);
        for k=j
            dec(:,:,i) = dec(:,:,i) & (O{i} > O{k}); 
        end  
    end
    
    img = zeros(1000,1000);
    
    decb = zeros(size(x1,2),size(x1,2));
    for i=1:size(dist.comp,2)
        decb = decb | dec(:,:,dist.comp(i));
    end
    
    img = 0.5 * (1 - decb) + (1-bwmorph(edge(decb,'canny'),'thicken'));
    
    subplot(1,2,2);imshow(img,[]);
    
    hold on; scatter(Y(:,1)*1000, Y(:,2)*1000, 10, 'x', 'b',...
                            'MarkerFaceColor',[0.2 0 1],...
                            'MarkerEdgeColor',[0.2 0 1]); 
                            
    grid on;
                        
    % Create xlabel
    xlabel({'','Principal Component - 1'});

    % Create ylabel
    ylabel({'Principal Component - 2',''});
    
    
    %%
    
    screen_size = get(0, 'ScreenSize');
    
    f = figure('units','normalized','outerposition',[0 0 1 1]);
        
    imshow(img,[]);
    set(gca, 'visible', 'on', 'XGrid', 'on', 'YGrid', 'on');

    hold on; scatter(Y(:,1)*1000, Y(:,2)*1000, 10, 'x', 'b',...
                            'MarkerFaceColor',[0.2 0 1],...
                            'MarkerEdgeColor',[0.2 0 1]);       
    
    
    title(fig_title);
    
    set(f, 'Position', [0 0 screen_size(3) screen_size(4) ] );
        
    tmp = getframe(gcf);
    imDB = tmp.cdata;
    
    close 
    %%

end


function [S] = classificationresult(dist, d, R)

    imH = d(1);
    imW = d(2);

    S = false(imH,imW);

    for i=1:size(dist.comp, 2)
        j = 1:size(R, 2);
        j = j(j ~= dist.comp(i));
        
        Sn = true(imH,imW);
        
        for k=1:size(j,2)
            Sn = Sn & (R{dist.comp(i)} > R{j(k)});
        end
        
        S = S | Sn;
    end
    
end


function [labels] = gmm_classify(data, obj, comp)
    
    [P] = posterior(obj, data);

    labels = false(size(data,1),1);
    
    for i=1:size(comp, 2)
        
        j = 1:size(P, 2);
        j = j(j ~= comp(i));
        
        Sn = true(size(data,1),1);
        
        for k=1:size(j,2)
            Sn = Sn & (P(:,comp(i)) > P(:,j(k)));
        end
        
        labels = labels | Sn;
    end
    
end


function [A] = decisionarea(obj, dist, d)
    
    imH = d(1);
    imW = d(2);
    
    x1 = 0.001:0.001:1;
    x2 = 0.001:0.001:1;

    [X1 X2] = meshgrid(x1,x2);

    [P] = posterior(obj, [X1(:) X2(:)]);

    for i=1:obj.NComponents
        R{i} = reshape(P(:,i), size(x1,2), size(x1,2));
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
    
    A = round(sum(decb(:) > 0)*1000 / (size(decb,1)*size(decb,2)))/10;
    
    % rp = regionprops(decb, 'area');
    % A = sum([rp.Area]);  
    
end


function [cva] = CrossValidation(Y, obj, comp)

    load fisheriris
    
    N = size(Y,1);
    indices = crossvalind('Kfold', N, 10);
    
    cp = classperf([ones(N-1,1);0]);
    
    for i = 1:10
        
        test = (indices == i); train = ~test;
        
        class = gmm_classify(Y(test, :), obj, comp);
        
        classperf(cp,class,test)
    end
    
    cva = cp.correctRate;
end
