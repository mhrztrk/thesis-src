%% LIBSVM Usage  
%
% 	"Usage: model = svmtrain(training_label_vector, training_instance_matrix, 'libsvm_options');\n"
% 	"libsvm_options:\n"
% 	"-s svm_type : set type of SVM (default 0)\n"
% 	"	0 -- C-SVC\n"
% 	"	1 -- nu-SVC\n"
% 	"	2 -- one-class SVM\n"
% 	"	3 -- epsilon-SVR\n"
% 	"	4 -- nu-SVR\n"
% 	"-t kernel_type : set type of kernel function (default 2)\n"
% 	"	0 -- linear: u'*v\n"
% 	"	1 -- polynomial: (gamma*u'*v + coef0)^degree\n"
% 	"	2 -- radial basis function: exp(-gamma*|u-v|^2)\n"
% 	"	3 -- sigmoid: tanh(gamma*u'*v + coef0)\n"
% 	"	4 -- precomputed kernel (kernel values in training_instance_matrix)\n"
% 	"-d degree : set degree in kernel function (default 3)\n"
% 	"-g gamma : set gamma in kernel function (default 1/num_features)\n"
% 	"-r coef0 : set coef0 in kernel function (default 0)\n"
% 	"-c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"
% 	"-n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)\n"
% 	"-p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)\n"
% 	"-m cachesize : set cache memory size in MB (default 100)\n"
% 	"-e epsilon : set tolerance of termination criterion (default 0.001)\n"
% 	"-h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)\n"
% 	"-b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)\n"
% 	"-wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)\n"
% 	"-v n : n-fold cross validation mode\n"
% 	"-q : quiet mode (no outputs)\n"
%

%% Workspace 

addpath('C:\Users\mahirztrk\Thesis\hasat\libsvm\bin',...
        'C:\Users\mahirztrk\Thesis\hasat\meanshift\bin',...
        'C:\Users\mahirztrk\Thesis\hasat\bilateral_filter\Bilateral Filtering',...
        'C:\Users\mahirztrk\Thesis\hasat\DrawLine',...
        'C:\Users\mahirztrk\Thesis\hasat\RoadTracking',...
        'C:\Users\mahirztrk\Thesis\hasat');


addpath('D:\TEZ_Calisma\hasat\libsvm\bin',...
        'D:\TEZ_Calisma\hasat\meanshift\bin',...
        'D:\TEZ_Calisma\hasat\bilateral_filter\Bilateral Filtering',...
        'D:\TEZ_Calisma\hasat\DrawLine',...
        'D:\TEZ_Calisma\hasat\RoadTracking',...
        'D:\TEZ_Calisma\hasat');    
    
%% SVM Training

road_vec = get_feature_set(image,road_mask);
nonroad_vec = get_feature_set(image,nonroad_mask);

Data = [road_vec;nonroad_vec];

rcnt = size(road_vec,1);
nrcnt = size(nonroad_vec,1);

groups = zeros(rcnt + nrcnt, 1);
groups(1:rcnt, :) = ones(rcnt, 1);
groups((rcnt+1):(rcnt + nrcnt), :) = -1 * ones(nrcnt, 1);
clear nrcnt cnt 

rmpath('C:\Program Files\MATLAB\R2010a\toolbox\bioinfo\biolearning');

% !!! Parameters to play with !!!
%   -g gamma    : set gamma in kernel function (default = 1/num_features)
%   -c cost     : set the parameter C (default = 1)
%   -wi weight  : set the parameter C of class i to weight*C (default = 1)
models = svmtrain(double(groups), double(Data), '-t 2 -g 8 -c 2'); 




%% Meanshift Segmentation 

%   'SpatialBandWidth' - segmentation spatial radius (integer) [7]
%   'RangeBandWidth'   - segmentation feature space radius (float) [6.5]
%   'MinimumRegionArea'- minimum segment area (integer) [20]
[fimage labels modes regSize grad conf]  = edison_wrapper(  image,... 
                                                            @RGB2Luv, ...
                                                            'SpatialBandWidth', 10, ...
                                                            'RangeBandWidth', 5, ...
                                                            'MinimumRegionArea', 150);

% colorvector = CalcColorVec(labels, image);
% segm_image = uint8(ColorizeSegments(labels, colorvector));
% or
modes = modes';

tmp = zeros(1, size(modes,1), size(modes,2));
tmp(1,:,:) = modes;
tmp = Luv2RGB(tmp);
modes = zeros(size(modes,1),size(modes,2));

modes(:,:) = tmp(1,:,:);
segm_image = Luv2RGB(fimage);


%figure;imagesc(segm_image);
figure;imagesc(labels);
figure;imshow(label2rgb(labels, @jet, 'k','shuffle'));

%% Segment based classification 

imH = size(image,1);
imW = size(image,2);
p = size(image,3);

rgn_cnt = max(max(labels)) + 1;

%%TestData = modes;
%image = zeros(1024,1024,5); image(:,:,1:4) = imagesc_5_1; image(:,:,5) = get_ndvi(image(:,:,1:4));

TestData = zeros(rgn_cnt,p);
for i=1:rgn_cnt
   [m n] = find(labels==(i-1), 1, 'first');
   TestData(i,:) = segm_image(m,n,:); 
end

test_label = ones(rgn_cnt,1);
[predicted_label, accuracy, decision_values] = svmpredict(double(test_label), double(256*TestData), models);

classified_image = zeros(imH, imW);

for i=1:rgn_cnt
    if(predicted_label(i) == 1)

        mask = (labels == i - 1);     
        classified_image = classified_image + double(255*mask);
    end
end

col_image = ColorizeClassifiedSegments(image, labels, predicted_label);
figure;imagesc(col_image);

figure;imagesc(classified_image);

%% Pixel based classification 

[imH imW p] = size(image);
TestData = reshape(image, imH*imW, p);
test_label = ones(size(TestData,1),1); 

[predicted_label, accuracy, decision_values] = svmpredict(double(test_label), double(TestData), models);
classified_image = reshape(predicted_label, imH, imW);

figure;imshow(classified_image,[]);

%% AdaBoost

[imH imW p] = size(image); TestData = reshape(image, imH*imW, p);                                                                  
[strong_h]=test_adaboost_single_set(TestData,feature_selection_matrix);
classified_image = reshape(strong_h, imH, imW);figure;imshow(classified_image,[]);   
 


%% Steerable Gauss Filter

% Compute directional filters
theta = 0:15:360;
for i = 1:length(theta)
   [J,H] = steerGauss([],theta(i),2,false);
   filters{i} = H;
end

% Using pre-computed filters.
J = -Inf*ones(size(I,1), size(I,2),2);
for i = 1:length(filters)
   [J(:,:,1),H] = steerGauss(I,filters{i},false);
   result = max(J,[],3); 
end

figure;imshow(result,[]);

%%
i = find(isnan(x) | isinf(x)); % Find bad elements in x
x(i) = []; % and delete them

% Alternatively,
i = find(~isnan(x) & ~isinf(x)); % Find elements that are not NaN and not infinite
x = x(i); % Keep those elements

% Both of these solutions can be further streamlined by using logical indexing:
x(isnan(x) | isinf(x)) = []; % Delete bad elements
% or
x = x(~isnan(x) & ~isinf(x)); % Keep good elements


% The sinc function has a piecewisede definition,
% This code uses find with vectorized computation to handle the two cases separately:

function y = sinc(x)
    % Computes the sinc function per element for a set of x values.
    y = ones(size(x)); % Set y to all ones, sinc(0) = 1
    i = find(x ~= 0); % Find nonzero x values
    y(i) = sin(x(i))./x(i); % Compute sinc where x ˜= 0
end

% A concise alternative is 
y = (sin(x) + (x==0)) ./ (x + (x==0));


min(min(A));

% A better method that works regardless of the number
[MinValue,MinIndex] = min(A(:)); % Find the minimum element in A
% The minimum value is MinValue, the index is MinIndex
MinSub = ind2sub(size(A),MinIndex); % Convert MinIndex to subscripts

