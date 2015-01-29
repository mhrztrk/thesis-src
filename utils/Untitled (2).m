for i=1:4
    mask{i} = bwmorph(ref{i}, 'skel', Inf);
end

img_name{1} = 'img1';
img_name{2} = 'img2';
img_name{3} = 'img3';
img_name{4} = 'img4';

for i=1:4
    model{i} = imgpca_2(uint8(img{i}*255), 'auto', 3);
end

for i=1:4
    fset{i,1}(:,:,1) = mat2gray(model{i}.scores(:,:,1),[0 255]) * 2 -1;
    fset{i,1}(:,:,2) = mat2gray(model{i}.scores(:,:,2),[0 255]) * 2 -1;
end


%%
for i=1:4
    
    for j=1:1
        nu = [0.01 0.05 0.1];
        range = [-1 1];

        for k=1:size(nu,2)
            [R] = imclassify(fset{i,j}, mask{i}, ref{i}, 'oc-svm', nu(k), range, img_name{i}, sprintf('Results_2018_01_16_ocsvm/%d',j));
            close all;
        end
    end
    
end

%%


for i=1:10
    [Rmin] = gmm_classification_ex((fset{i,1}*0.5 + 0.5), mask{i}, ref{i}, img_name{i}, 'Results_2012_12_13/gmm4');
    close all
end

%%
for i=1:4
    clres_eval{i}(:,:,1) = double(            (ref{i}>0) & imcomplement(clres_ocsvm{i}));
    clres_eval{i}(:,:,2) = double(            (ref{i}>0) &             (clres_ocsvm{i}));
    clres_eval{i}(:,:,3) = double(imcomplement(ref{i}>0) &             (clres_ocsvm{i}));
    
    figure;imshow(clres_eval{i},[])
    
    imwrite(clres_eval{i},sprintf('ocsvm_clres_eval_img_%d.png',i));
end
%%