%%
mask{ 1} = bwmorph(im_14_mask, 'skel', Inf);
mask{ 2} = bwmorph(im_17_mask, 'skel', Inf);
mask{ 3} = bwmorph(im_27_mask, 'skel', Inf);
mask{ 4} = bwmorph(im_4_mask, 'skel', Inf);
mask{ 5} = bwmorph(im_2_1_mask, 'skel', Inf);
mask{ 6} = bwmorph(im_19_mask, 'skel', Inf);
mask{ 7} = bwmorph(img1_mask, 'skel', Inf);
mask{ 8} = bwmorph(img2_mask, 'skel', Inf);
mask{ 9} = bwmorph(img3_mask, 'skel', Inf);
mask{10} = bwmorph(img4_mask, 'skel', Inf);

ref{ 1} = im_14_mask; 
ref{ 2} = im_17_mask;
ref{ 3} = im_27_mask;
ref{ 4} = im_4_mask;
ref{ 5} = im_2_1_mask;
ref{ 6} = im_19_mask;
ref{ 7} = img1_mask;
ref{ 8} = img2_mask;
ref{ 9} = img3_mask;
ref{10} = img4_mask;

img_name{ 1} = 'im14';
img_name{ 2} = 'im17';
img_name{ 3} = 'im27';
img_name{ 4} = 'im4';
img_name{ 5} = 'im2_1';
img_name{ 6} = 'im19';
img_name{ 7} = 'img1';
img_name{ 8} = 'img2';
img_name{ 9} = 'img3';
img_name{10} = 'img4';

fset{ 1,1}(:,:,1) = mat2gray(model14.scores(:,:,1),[0 255]) * 2 -1;
fset{ 1,1}(:,:,2) = mat2gray(model14.scores(:,:,2),[0 255]) * 2 -1;

fset{ 1,2}(:,:,1) = im14_pfilt(:,:,1) * 2 - 1;
fset{ 1,2}(:,:,2) = im14_pfilt(:,:,2) * 2 - 1;
fset{ 1,2}(:,:,3) = im14_pfilt(:,:,3) * 2 - 1;
fset{ 1,2}(:,:,4) = im14_pfilt(:,:,4) * 2 - 1;

fset{ 2,1}(:,:,1) = mat2gray(model17.scores(:,:,1),[0 255]) * 2 -1;
fset{ 2,1}(:,:,2) = mat2gray(model17.scores(:,:,2),[0 255]) * 2 -1;

fset{ 2,2}(:,:,1) = im17_pfilt(:,:,1) * 2 - 1;
fset{ 2,2}(:,:,2) = im17_pfilt(:,:,2) * 2 - 1;
fset{ 2,2}(:,:,3) = im17_pfilt(:,:,3) * 2 - 1;
fset{ 2,2}(:,:,4) = im17_pfilt(:,:,4) * 2 - 1;

fset{ 3,1}(:,:,1) = mat2gray(model27.scores(:,:,1),[0 255]) * 2 -1;
fset{ 3,1}(:,:,2) = mat2gray(model27.scores(:,:,2),[0 255]) * 2 -1;

fset{ 3,2}(:,:,1) = im27_pfilt(:,:,1) * 2 - 1;
fset{ 3,2}(:,:,2) = im27_pfilt(:,:,2) * 2 - 1;
fset{ 3,2}(:,:,3) = im27_pfilt(:,:,3) * 2 - 1;
fset{ 3,2}(:,:,4) = im27_pfilt(:,:,4) * 2 - 1;

fset{ 4,1}(:,:,1) = mat2gray(model4.scores(:,:,1),[0 255]) * 2 -1;
fset{ 4,1}(:,:,2) = mat2gray(model4.scores(:,:,2),[0 255]) * 2 -1;

fset{ 4,2}(:,:,1) = im4_pfilt(:,:,1) * 2 - 1;
fset{ 4,2}(:,:,2) = im4_pfilt(:,:,2) * 2 - 1;
fset{ 4,2}(:,:,3) = im4_pfilt(:,:,3) * 2 - 1;
fset{ 4,2}(:,:,4) = im4_pfilt(:,:,4) * 2 - 1;
 
fset{ 5,1}(:,:,1) = mat2gray(model2_1.scores(:,:,1),[0 255]) * 2 -1;
fset{ 5,1}(:,:,2) = mat2gray(model2_1.scores(:,:,2),[0 255]) * 2 -1;

fset{ 5,2}(:,:,1) = im_2_1(:,:,1) * 2 - 1;
fset{ 5,2}(:,:,2) = im_2_1(:,:,2) * 2 - 1;
fset{ 5,2}(:,:,3) = im_2_1(:,:,3) * 2 - 1;
fset{ 5,2}(:,:,4) = im_2_1(:,:,4) * 2 - 1;

fset{ 6,1}(:,:,1) = mat2gray(model19.scores(:,:,1),[0 255]) * 2 -1;
fset{ 6,1}(:,:,2) = mat2gray(model19.scores(:,:,2),[0 255]) * 2 -1;

fset{ 6,2}(:,:,1) = im_19(:,:,1) * 2 - 1;
fset{ 6,2}(:,:,2) = im_19(:,:,2) * 2 - 1;
fset{ 6,2}(:,:,3) = im_19(:,:,3) * 2 - 1;
fset{ 6,2}(:,:,4) = im_19(:,:,4) * 2 - 1;

fset{ 7,1}(:,:,1) = mat2gray(model_img1.scores(:,:,1),[0 255]) * 2 -1;
fset{ 7,1}(:,:,2) = mat2gray(model_img1.scores(:,:,2),[0 255]) * 2 -1;

fset{ 7,2}(:,:,:) = img1 * 2 - 1;

fset{ 8,1}(:,:,1) = mat2gray(model_img2.scores(:,:,1),[0 255]) * 2 -1;
fset{ 8,1}(:,:,2) = mat2gray(model_img2.scores(:,:,2),[0 255]) * 2 -1;

fset{ 8,2}(:,:,:) = img2 * 2 - 1;

fset{ 9,1}(:,:,1) = mat2gray(model_img3.scores(:,:,1),[0 255]) * 2 -1;
fset{ 9,1}(:,:,2) = mat2gray(model_img3.scores(:,:,2),[0 255]) * 2 -1;

fset{ 9,2}(:,:,:) = img3 * 2 - 1;

fset{10,1}(:,:,1) = mat2gray(model_img4.scores(:,:,1),[0 255]) * 2 -1;
fset{10,1}(:,:,2) = mat2gray(model_img4.scores(:,:,2),[0 255]) * 2 -1;

fset{10,2}(:,:,:) = img4 * 2 - 1;

%%
for i=7:10
    
    for j=1:2
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