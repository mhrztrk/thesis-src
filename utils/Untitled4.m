
x1 = 0.001:0.001:1;
x2 = 0.001:0.001:1;

[X1 X2] = meshgrid(x1,x2);

F1 = mvnpdf([X1(:) X2(:)],obj.mu(1,:),obj.Sigma(:,:,1));
F2 = mvnpdf([X1(:) X2(:)],obj.mu(2,:),obj.Sigma(:,:,2));

F = reshape(F1,length(x2),length(x1)) + reshape(F2,length(x2),length(x1));
surf(x1,x2,F);

O1 = pdf(gmdistribution(obj.mu(1,:), obj.Sigma(:,:,1)),[X1(:) X2(:)]);
O2 = pdf(gmdistribution(obj.mu(2,:), obj.Sigma(:,:,2)),[X1(:) X2(:)]);

figure;scatter(Y(:,1), Y(:,2), 10, 'o', 'b',...
     'MarkerFaceColor',[0.2 0 1],...
    'MarkerEdgeColor',[0.2 0 1]); 

hold on
contour(x1,x1,reshape(O1,size(x1,2),size(x1,2)));
contour(x1,x1,reshape(O2,size(x1,2),size(x1,2)));

% Create xlabel
xlabel({'','Principal Component - 1'});

% Create ylabel
ylabel({'Principal Component - 2',''});