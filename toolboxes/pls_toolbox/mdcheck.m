function [flag,missmap,data] = mdcheck(data,options);
%MDCHECK Missing Data Checker and infiller.
%  This function checks for missing data and infills it using a PCA model
%  if desired. The input is the data to be checked (data) as either a
%  double array or a dataset object. Optional input (options) is a
%  structure containing options for how the function is to run.
%  Outputs are the fraction of missing data (flag), a map of the locations
%  of the missing data as an unint8 variable (missmap), and the data with
%  the missing values filled in (infilled). Depending on the plots option,
%  a plot of the missing data may also be output.
%  
%  INPUTS:
%     data    : data to be checked.
%     options : optional structure of settings for function. 
%
%  OUTPUTS:
%     flag    : fraction of missing data.
%     mismap  : map of missing data (unint8).
%     data    : data with missing values filled in. 
%
%  OPTIONS:
%
%  options = stucture array with the following fields:
%          plots: [{'none'}| 'final' ] governs plot of missing data map.
%                 Map contains green pixels where data is OK, Yellow where
%                 data is missing but can be replaced, and Red where too
%                 much is missing to be replaced.
%       frac_ssq: [{0.95}] desired fraction (between 0 and 1) of variance
%                  to be captured by the PCA model,
%        max_pcs: [{5}] maximum number of PCs in model, if 0, then it uses
%                  the mean, 
%     meancenter: ['no' | {'yes'}], tells whether to use mean centering in
%                 the algorithm,
%     recalcmean: [{'no'} | 'yes'], recalculate mean centering after each
%                  cycle of replacement (may improve results for small
%                  matricies)
%        display: [{'off'} | 'on'] governs level of display,
%      tolerance: [{1e-6  100}] convergence criteria, the first element is
%                  the minimum change and the second is the maximum iterations,
%    max_missing: [{.4}] maximum fraction of missing data with which mdcheck 
%                  will operate, and
%        toomuch: [{'error'}| 'exclude' ] What action should be taken if too much 
%                  missing data is found. 'error' exit with error message, 
%                  'exclude' will exclude elements (rows/columns/slabs/etc)
%                  which contain too much missing data from the data before
%                  replacement. 'exclude' requires a dataset object as
%                  input for (data).
%      algorithm: [ {'svd'} | 'nipals' ] specifies missing data replacement
%                  algorithm to use. NIPALS typically used for large amounts of 
%                  missing data or large multi-way arrays.
%
%Note: MDCHECK captures up to options.frac_ssq of the variance using options.max_pcs
%      or fewer PCA components.
%
%I/O: [flag,missmap,infilled] = mdcheck(data,options);
%
%See also: PARAFAC, PCA

%Copyright Eigenvector Research, Inc. 2002-2008
%Licensee shall not re-compile, translate or convert "M-files" contained
% in PLS_Toolbox for use with any software other than MATLAB®, without
% written permission from Eigenvector Research, Inc.
%by BMW 5/9/02, debugged 5/16/02 BMW
%7/22/02 RB exchanged the unfolding of nway arrays so that as few data are missing per column as possible
%7/22/02 RB removed replace.m with a single imputation step which is faster
%7/22/02 RB removed warning for max iterations unless disp is on
%8/15/02 RB fixed 2D matrix error
%8/15/02 RB exchanged svd with pca
%26/10/02 RB exchanged PCA with NIPALS
%2/20/04 JMS reconcile options with reconopts for EVRIIO compatibility
%2/24/04 JMS added recalcmean option to recalculate the mean centering during replacement
%04/27/05 RSK Add I/O outline to help.
%26/7/05 RB modified NIPALS to catch division by zero for extreme data


if nargin == 0; data = 'io'; end
if ischar(data);
  options.name    = 'options';
  options.plots   = 'none';
  options.max_pcs = 5;
  options.frac_ssq = 0.95;
  options.meancenter = 'yes';
  options.recalcmean = 'no';
  options.display = 'off';
  options.tolerance = [1e-6 100];
  options.max_missing = .4;
  options.toomuch = 'error'; % 'error' | 'exclude'
  options.algorithm = 'svd';
  
  if nargout==0 
    clear flag; 
    evriio(mfilename,data,options); 
  else 
    flag = evriio(mfilename,data,options); 
  end
  return 
end

if nargin < 2;
  options = [];
end
options = reconopts(options,'mdcheck');

%if asked to, exclude elements with too much missing data
if strcmp(options.toomuch,'exclude')
  data = excludemissing(data,options.max_missing);
end

%if is dataset, explode into non-dso
wasdso = isa(data,'dataset');
if wasdso
  origdata = data;
  incl = data.include;
  data = data.data(incl{:});
end

flag = mean(data); 
while length(flag)>1; flag = mean(flag); end; 
flag = ~isfinite(flag);

if ~flag
  missmap = [];
else
  missmap = uint8(~isfinite(data));
  
  flag    = sum(missmap);
  while length(flag)>1; flag = sum(flag); end;
  flag    = flag./prod(size(data));           %calculate fraction of missing data
  
  if nargout > 2
    datasize = size(data);
    datadims = ndims(data);
    if datadims == 2 & min(datasize)==1;
      disp('Can not infill vector data')
      infilled = [];
    elseif datadims >= 2
      
      if datadims>2
        % RB added. Try other unfoldings to find the one with lowest percentage missing in any row or column
        id = [2:datadims 1];
        for dd=1:datadims;
          data2 = permute(data,[dd:datadims 1:dd-1]);
          datasize2 = size(data2);
          data2 = reshape(data2,prod(datasize2(1:floor(datadims/2))),prod(datasize2(ceil(datadims/2+.1):end)));
          missmap = uint8(~isfinite(data2));
          lowmean(dd) = max(max(mean(double(missmap'))),max(mean(double(missmap))));
        end
        [alow,blow]=min(lowmean);
        clear data2
        data = permute(data,[blow(1):datadims 1:blow(1)-1]);
        datasize = size(data);
        data = reshape(data,prod(datasize(1:floor(datadims/2))),prod(datasize(ceil(datadims/2+.1):end)));
      else
        datasize = size(data);
      end
      missmap = uint8(~isfinite(data));
      [m,n] = size(data);
      % Check for any columns or rows with too much missing
      rowmean = mean(double(missmap'));
      colmean = mean(double(missmap));
      if any(rowmean>options.max_missing) | any(colmean>options.max_missing)
        error('Too much missing data in some rows or columns of unfolded data')
      end
      % Start by replacing with mean of available data
      for i = 1:n
        use = isfinite(data(:,i));
        if any(use);
          mx = mean(data(use,i));
        else
          mx = 0;
        end
        data(~isfinite(data(:,i)),i) = mx;
      end
      if options.max_pcs > 0
        switch lower(options.meancenter)
          case {'yes'}
            [data,mx] = mncn(data);
          case {'no'}
            % Do nothing
          otherwise 
            disp('unknown centering option, no centering used')
        end  
        % Replace based on options PCs
        totssq = sum(sum(data.^2));
        pcs = 0;
        frac_ssq = 0;
        while pcs < options.max_pcs & frac_ssq < options.frac_ssq
          pcs = pcs + 1;
          change = 1;
          count  = 0;
          while change > options.tolerance(1)
            count = count + 1;
            if strcmp(lower(options.display),'on')  
              disp(sprintf('Now working on iteration number %g',count));
            end
            
            switch lower(options.algorithm)
              case 'nipals'
                x = data;
                x(find(missmap)) = NaN; % Put back NaNs (should be changed to not infilling them at all when nipals is chosen)
                [scores,loads] = pcanipals(x,pcs,1);
                %               case 'pca'
                %                 [ssq,datarank,loads,scores,msg] = pcaengine(data,pcs,pcaopt);
              case 'svd'
                if n < m
                  cov = (data'*data)/(m-1);
                  [u,s,v] = svd(cov);
                else
                  cov = (data*data')/(m-1);
                  [u,s,v] = svd(cov);  
                  v = data'*v;
                  for i = 1:pcs
                    v(:,i) = v(:,i)/norm(v(:,i));
                  end
                end
                loads = v(:,1:pcs);
                scores = data*loads;
                
              otherwise
                error('Unrecognized algorithm')
            end
            
            estdata = data;
            mddd=scores*loads';
            estdata(find(missmap)) = mddd(find(missmap));
            
            if strcmp(lower(options.recalcmean),'yes') & strcmp(lower(options.meancenter),'yes')
              %Re-evaluate mean-centering for new data estimation
              [estdata,mx] = mncn(rescale(estdata,mx));     %re-evaluate mean-centering
            end  
            
            dif = data - estdata;
            change = sum(sum(dif.^2));
            if strcmp(lower(options.display),'on')
              disp(sprintf('Sum of squared differences in missing data estimates = %g',change));
            end
            data = estdata;
            if count == options.tolerance(2);
              if strcmp(lower(options.display),'on')  
                disp(['Algorithm failed to converge after ',int2str(options.tolerance(2)),' iterations.'])
              end
              change = 0;
            end
          end
          if strcmp(lower(options.display),'on')
            disp('  ')
            disp('Now forming final PCA model')
            disp('   ')
          end
          ssq = [[1:pcs]' zeros(pcs,2)];
          for i = 1:pcs
            resmat = data - scores(:,1:i)*loads(:,1:i)';
            resmat = resmat - resmat.*double(missmap>0);
            ssqres = sum(sum(resmat.^2));
            ssq(i,3) = (1 - ssqres/totssq)*100;
          end
          ssq(1,2) = ssq(1,3);
          for i = 2:pcs
            ssq(i,2) = ssq(i,3) - ssq(i-1,3);
          end
          if strcmp(lower(options.display),'on')
            disp('    Percent Variance Captured') 
            disp('           by PCA Model')
            disp('     Based on Known Data Only')
            disp('  ')
            disp('    PC#      %Var      %TotVar')
            disp(ssq)
          end  
          frac_ssq = ssq(end,3)/100;
        end
      end
      switch lower(options.meancenter)
        case {'yes'}
          data = rescale(data,mx);
        case {'no'}
          % Do nothing
      end 
      data(missmap==2) = NaN;  %replace "too much missing" with NaN's again
      if datadims > 2%reshape back to multiway data
        data = reshape(data,datasize);
        data = ipermute(data,[blow(1):datadims 1:blow(1)-1]);
        missmap = reshape(missmap,datasize);
        missmap = ipermute(missmap,[blow(1):datadims 1:blow(1)-1]);
      end
    end
    
    %if was DSO, adjust missmap size and replaced data for original data
    if wasdso
      missmap = nassign(zeros(size(origdata)),missmap,incl,1:length(incl));
      origdata.data = nassign(zeros(size(origdata)),data,incl,1:length(incl));
      data = origdata;  %pass back replaced
    end

  end
end

if strcmp(options.plots,'final')
  immap = single(missmap);
  if isempty(immap)
    immap = zeros(size(data));
  end
  %if was a DSO, hide excluded elements in image
  if wasdso
    for mode = 1:length(incl)
      excl = setdiff(1:size(immap,mode),incl{mode});
      immap = nassign(immap,nan,excl,mode);
    end
  end

  %unfold map into two-way
  immap = immap(:,:);
  %and classify elements as not missing, OK missing, or Too Much missing
  for mode = 1:2;
    bad  = find(mean(immap,3-mode)>options.max_missing);
    badx = nindex(immap,bad,mode);
    badx(badx==1) = 2;
    immap = nassign(immap,badx,bad,mode);
  end
  
  immap(isnan(immap)) = -1;  %set excluded items to -1

  %do plot with green, yellow, red coloration for missing data
  fig = figure;
  colormap([1 1 1; 0 .4 0; 1 1 0; 1 0 0])
  imagesc(immap);
  caxis([-1 2]);
  title('Missing Data Map (Show legend for key)')
  ylabel('Row');
  if ndims(data) > 2
    xlabel('Other Modes Unfolded');
  else
    xlabel('Column');
  end
  
  hold on
  h = plot(nan,1,'ys',nan,1,'rs');
  set(h(1),'markerfacecolor','yellow')
  set(h(2),'markerfacecolor','red')
  legendname(h(1),'Missing - Can Replace');
  legendname(h(2),'Missing - Can Not Replace');
  hold off
  
  set(fig,...
    'windowbuttonmotionfcn','showposition(gcbf,1)',...
    'pointer','crosshair');
end


%--------------------------------------------------------
function [t,p,Mean,Fit,RelFit] = pcanipals(X,F,cent);

% NIPALS-PCA WITH MISSING ELEMENTS
% 20-6-1999
%
% Calculates a NIPALS PCA model. Missing elements 
% are denoted NaN. The solution is nested
% 
% Comparison for data with missing elements
% NIPALS : Nested    , not least squares, not orthogonal solutoin
% LSPCA  : Non nested, least squares    , orthogonal solution
% 
% I/O
% [t,p,Mean,Fit,RelFit] = pcanipals(X,F,cent);
% 
% X   : Data with missing elements set to NaN
% F   : Number of componets
% cent: One if centering is to be included, else zero
% 
% Copyright
% Rasmus Bro
% KVL 1999
% rb@kvl.dk

% Oct, 2002, rb, modified to allow smaller number of iterations only and to have
% bi display
[I,J]=size(X);
if any(sum(isnan(X))==I)|any(sum(isnan(X)')==J)
  error(' One column or row only contains missing')
end

Xorig      = X;
Miss       = isnan(X);
NotMiss    = ~isnan(X);
ssX    = sum(X(find(NotMiss)).^2);

Mean   = zeros(1,J);
if cent
  Mean    = nanmean(X);
end
X      = X - ones(I,1)*Mean;

t=[];
p=[];

for f=1:F
  Fit    = 3;
  OldFit = 6;
  it     = 0;
  
  T      = nanmean(X')';
  P      = nanmean(X)';
  if std(T)<eps*1000;
       T = T + randn(size(T));
   end
   if std(P)<eps*1000;
       P = P + randn(size(P));
   end

  Fit    = 2;
  FitOld = 3;
  
  while abs(Fit-FitOld)/FitOld>1e-7 & it < 50;
    FitOld  = Fit;
    if FitOld ==0;
      FitOld = eps;
    end
    it      = it +1;
    for j = 1:J
      id=find(NotMiss(:,j));
      if length(id)>0
        if T(id)'*T(id)~=0
          P(j) = T(id)'*X(id,j)/(T(id)'*T(id));
        else
          P(j)=0;
        end
      else
        P(j)=0;
      end
    end
    normP=norm(P);
    if normP~=0
      P = P/norm(P);
    else
      P = P;
    end
    for i = 1:I
      id=find(NotMiss(i,:));
      if length(id)>0
        if P(id)'*P(id)~=0
          T(i) = P(id)'*X(i,id)'/(P(id)'*P(id));
        else
          T(i)=0;
        end
      else
        T(i)=0;
      end
    end
    Fit = X-T*P';
    Fit = sum(Fit(find(NotMiss)).^2);
  end
  t = [t T];
  p = [p P];
  X = X - T*P';
end

Model   = t*p' + ones(I,1)*Mean;
Fit     = sum(sum( (Xorig(find(NotMiss)) - Model(find(NotMiss))).^2));
RelFit  = 100*(1-Fit/ssX);

function y = nanmean(x)
if isempty(x) % Check for empty input.
  y = NaN;
  return
end

% Replace NaNs with zeros.
nans = isnan(x);
i = find(nans);
x(i) = zeros(size(i));

if min(size(x))==1,
  count = length(x)-sum(nans);
else
  count = size(x,1)-sum(nans);
end

% Protect against a column of all NaNs
i = find(count==0);
count(i) = ones(size(i));
y = sum(x)./count;
y(i) = i + NaN;
