function varargout = pca(varargin)
%PCA Principal components analysis.
%  PCA uses sigular value decomposition on a data matrix (X) and
%  returns scores (T) and loadings (P) which describe the data matrix
%     X = TP'+E
%
%  INPUTS:
%         x = X-block (2-way array class "double" or "dataset"), and
%     ncomp = number of components to to be calculated (positive integer scalar).
%
%  OPTIONAL INPUT:
%   options = structure array with the following fields:
%           display: [ 'off' | {'on'} ]      governs level of display to command window.
%             plots: [ 'none' | {'final'} ]  governs level of plotting.
%     outputversion: [ 2 | {3} ]             governs output format.
%         algorithm: [ {'svd'} | 'maf' | 'robustpca' ]     algorithm for decomposition.
%                      Algorithm 'maf' requires Eigenvector's MIA_Toolbox.
%     preprocessing: { [] }                  preprocessing structure (see PREPROCESS).
%      blockdetails: [ 'compact' | {'standard'} | 'all' ]   Extent of predictions and raw residuals  
%                     included in model. 'standard' = none, 'all' x-block 
%   confidencelimit: [{0.95}] Confidence level for Q and T2 limits. A value
%                     of zero (0) disables calculation of confidence
%                     limits.
%          roptions: structure of options to pass to robpca (robust PCA engine from the Libra
%                    Toolbox). 
%                alpha: (1-alpha) measures the number of outliers the
%                       algorithm should resist. Any value between 0.5 and 1 may
%                       be specified (default = 0.75). These options are
%                       only used when algorithm is robust PCA.
%
%  OUTPUT:
%     model = standard model structure containing the PCA model (See MODELSTRUCT)
%
%I/O: model   = pca(x,ncomp,options);  %identifies model (calibration step)
%I/O: pred    = pca(x,model,options);  %projects a new X-block onto existing model
%I/O: options = pca('options');        %returns default options structure
%I/O: pca demo                         %runs a demo of the PCA function.
%
%See also: ANALYSIS, EVOLVFA, EWFA, EXPLODE, PARAFAC, PLOTLOADS, PLOTSCORES, PREPROCESS, SSQTABLE

%Copyright Eigenvector Research, Inc. 1991-2008
%Licensee shall not re-compile, translate or convert "M-files" contained
% in PLS_Toolbox for use with any software other than MATLAB®, without
% written permission from Eigenvector Research, Inc.
%JMS 3/7/02 -modified PLS shell for use in PCA
%jms 3/19/02 -added more support for model predictions
%jms 3/22/02 -fixed variance captured 100% bug
%jms 8/6/03 - fixed extra "see also" help

%Old I/O:
%  The input is the data matrix (data). Outputs are the scores
%  (scores), loadings (loads), variance info (ssq), residuals
%  (res), Q limit (reslm), T^2 limit (tsqlm), and T^2's (tsq).
%
%  Optional inputs are (plots) plots = 0 suppresses all plots,
%  plots =  1 [default] produces plots with no confidence limits,
%  plots =  2 produces plots with limits,
%  plots = -1 plots the eigenvalues only (without limits),
%  a vector (scl) for plotting scores against, (if scl = 0 sample 
%  numbers will be used), and a scalar (lv) which specifies the
%  number of principal components to use in the model and
%  which suppresses the prompt for number of PCs.
%
%I/O: [scores,loads,ssq,res,reslm,tsqlm,tsq] = pca(data,plots,scl,lvs);
%
%  Note: with plots = 0 and lv specified, this routine requires
%  no interactive user input. If you would like to scale the data
%  before processing use the functions AUTO or SCALE.
%

%Start Input
if nargin==0  % LAUNCH GUI
  analysis pca
  return
elseif ischar(varargin{1}) %Help, Demo, Options

  options = [];
  options.name          = 'options';
  options.display       = 'on';     %Displays output to the command window
  options.plots         = 'final';  %Governs plots to make
  options.outputversion = 3;        %2,3 Tells what to output (3=ModelStruct)
  options.preprocessing = {[]};     %See Preprocess
  options.algorithm     = 'svd';    %SVD algorithm
  options.blockdetails  = 'standard';  %level of details
  options.confidencelimit = 0.95;
  options.roptions.alpha = 0.75; %Alpha for robust methods between 0.5 and 1.
  options.definitions   = @optiondefs;

  if nargout==0; evriio(mfilename,varargin{1},options); else; varargout{1} = evriio(mfilename,varargin{1},options); end
  return; 
  
else
  if nargin<2
    error([ upper(mfilename) ' requires 2 inputs.'])
  end
  
  %A) Check Options Input
  predictmode = 0;    %default is calibrate mode

  %options if we interpret as an old call
  v2options.display       = 'on';         %Displays output to the command window
  v2options.plots         = 'final';      %Governs plots to make
  v2options.outputversion = 2;            %2,3 Tells what to output (2=[b,ssq,t,p,eigs])
  v2options.preprocessing = {[]};         %See preprocess
  %v2options.algorithm     = 'svd';        %SVD algorithm
  v2options.blockdetails  = 'standard';  %level of details
  
  switch nargin
  case 2  %two inputs
    %v3 : (x,model)
    %v3 : (x,ncomp)
    %v2 : (x,plots)  (invalid - used to be valid in pre v3.0)

    if isa(varargin{2},'struct')
      %v3 : (x,model)
      varargin = {varargin{1}, [], pca('options'), varargin{2}};
    elseif varargin{2}>0;
      %v3 : (x,ncomp)
      varargin{3} = pca('options');
    else
      %v2 : (x,plots)
      error(['Input NCOMP is missing. Type: ''help ' mfilename '''']);      
    end
        
  case 3  %three inputs
    %v3 : (x,model,options)
    %v3 : (x,ncomp,options)    
    %v2 : (x,plots,scl)      (invalid - used to be valid in pre v3.0)
    
    if isa(varargin{3},'struct');
      if isa(varargin{2},'struct');
        %v3 : (x,model,options)
        varargin = {varargin{1}, [], varargin{3}, varargin{2}};
      else
        %v3 : (x,ncomp,options)    
      end
    else
      %v2 : (x,plots,scl)      (invalid - used to be valid in pre v3.0)
      error(['Input NCOMP is missing. Type: ''help ' mfilename ''''])
    end
    
  case 4   %four inputs
    %v2 : (x,plots,scl,ncomp)
    %v3 : (x,ncomp,model,options)

    if ~isa(varargin{4},'struct');
      %v2 : (x,plots,scl,ncomp)
      switch varargin{2}
      case 0
        v2options.display       = 'off';         %Displays output to the command window
        v2options.plots         = 'none';        %Governs plots to make
      otherwise
        error(['Input OPTIONS or MODEL not recognized. Type: ''help ' mfilename ''''])
      end
      varargin = {varargin{1} varargin{4} v2options};
      
    elseif isa(varargin{3},'struct');
      %v3 : (x,ncomp,model,options)
      %v3 : (x,ncomp,options,model)   (technically invalid but we'll accept it anyway)
      if ~isfield(varargin{4},'modeltype');
        varargin([3 4]) = varargin([4 3]);
      end
      
    else
      error(['Input MODEL not recognized. Type: ''help ' mfilename ''''])
    end

  otherwise
    error(['Unrecognized input. Type: ''help ' mfilename ''''])        
  end

  try
    options = reconopts(varargin{3},pca('options'));
  catch
    error(['Input OPTIONS not recognized. Type: ''help ' mfilename ''''])
  end
  
  if ~ismember(options.algorithm,{'svd' 'maf' 'robustpca'})
    error(['Algorithm [' options.algorithm '] not recognized. PCA supports ''svd'', ''maf'', and ''robustpca''.'])
    return
  end
  
  switch options.outputversion
  case{2,3}
    %Take no action these are ok inputs
  otherwise
    options.outputversion = 3;
    warning('OPTIONS.OUTPUTVERSION not recognized. Reset to 3.')
  end
  options.blockdetails = lower(options.blockdetails);
  if ~ismember(options.blockdetails,{'compact','standard','all'})
    error(['OPTIONS.BLOCKDETAILS not recognized. Type: ''' mfilename ' options'''])
  end
  if ~isfield(options,'rawmodel'); 
    options.rawmodel = 0;       %undocumented option to output raw results ONLY
  end
  
  %B) check model format
  if length(varargin)>=4;
    try
      varargin{4} = updatemod(varargin{4});        %make sure it's v3.0 model
    catch
      error(['Input MODEL not recognized. Type: ''help ' mfilename ''''])
    end
    predictmode = 1;                                  %and set predict mode flag
    if isempty(varargin{2});
      varargin{2} = size(varargin{4}.loads{2,1},2);   %get ncomp from model (if needed)
    end
  end
  
  %C) CHECK Data Inputs
  datasource = {getdatasource(varargin{1})};
  if isa(varargin{1},'double')      %convert varargin{1} and varargin{2} to DataSets
    varargin{1}        = dataset(varargin{1});
    varargin{1}.name   = inputname(1);
    varargin{1}.author = 'PCA';
  elseif ~isa(varargin{1},'dataset')
    error(['Input X must be class ''double'' or ''dataset''.'])
  end
  if ndims(varargin{1}.data)>2
    error(['Input X must contain a 2-way array. Input has ',int2str(ndims(varargin{1}.data)),' modes.'])
  end
  
  %D) Check Meta-Parameters Input
  if isempty(varargin{2}) | prod(size(varargin{2}))>1 | varargin{2}<1 | varargin{2}~=fix(varargin{2});
    error('Input NCOMP must be integer scalar.')
  end
  if ~options.rawmodel & predictmode & varargin{2}~=size(varargin{4}.loads{2,1},2);
    error('Cannot use a different number of components (NCOMP) with previously created model');
  end

  %---------------------------------------------------------------
  % ready to go...
  x = varargin{1};
  
  %Handle Preprocessing
  if isempty(options.preprocessing);
    options.preprocessing = {[]};  %reinterpet as empty cell
  end
  if ~isa(options.preprocessing,'cell');
    options.preprocessing = {options.preprocessing};  %insert into cell
  end

  preprocessing = options.preprocessing;

  %Call PCA Function
  if ~predictmode | options.outputversion ~= 3;
    %basic calibrate call

    %Initialize model structure
    model      = modelstruct('pca');
    model.date = date;
    model.time = clock;

    model = copydsfields(x,model,[],{1 1});      %copy all mode labels, etc.
    model.datasource = datasource;
    
    if mdcheck(x.data(x.includ{1},x.includ{2}));
      if strcmp(options.display,'on'); warning('Missing Data Found - Replacing with "best guess". Results may be affected by this action.'); end
      [flag,missmap,replaced] = mdcheck(x.data(x.includ{1},x.includ{2}));
      x.data(x.includ{1},x.includ{2}) = replaced;
    end
    
    %preprocessing
    if ~isempty(preprocessing{1});
      [xpp,preprocessing{1}] = preprocess('calibrate',preprocessing{1},x);
    else
      xpp = x;
    end

    switch options.algorithm
      case {'svd' 'maf'}
        opts         = pcaengine('options');
        opts.display = 'off';
        [model.detail.ssq,datarank,model.loads{2,1}] = pcaengine(xpp.data(xpp.includ{1},xpp.includ{2}),varargin{2},opts);

        if strcmp(lower(options.algorithm),'maf')
          if exist('maf')
            [junk,model.loads{2,1},mn,model.detail.ssq] = maf(xpp.data(xpp.includ{1},xpp.includ{2}));
            model.loads{2,1} = model.loads{2,1}(:,1:min(end,min(datarank,varargin{2})));
          else
            error('Algorithm "maf" requires MIA_Toolbox');
          end
        end
      case 'robustpca'
        %Robust PCA

        if varargin{2}>50;
          error('Robust PCA can not be used for more than 50 principal components.');
        end

        robpcaout = robpca(xpp.data(xpp.includ{1},xpp.includ{2}),'k',varargin{2},'kmax',varargin{2},'alpha',options.roptions.alpha,'plots',0);
        model.detail.ssq = [[1:length(robpcaout.L)]' robpcaout.L(:) zeros(length(robpcaout.L),1) zeros(length(robpcaout.L),1)];
        datarank = robpcaout.k; %Number of (chosen) principal components
        model.loads{2,1} = robpcaout.P; %Robust loadings (eigenvectors)
        model.detail.robustpca = robpcaout;

        %Scale step, correct for Robust centering of model. Add to options and model.detail.prepro.
        xpp.data(:,xpp.includ{2}) = scale(xpp.data(:,xpp.includ{2}), model.detail.robustpca.M);
        %Add a prepro step to account for what Robust method does.
        robpp = preprocess('default','mean center');
        robpp.description = 'Robust Centering';
        robpp.tooltip = 'Remove robustly calculated offset from each variable.';
        robpp.keyword = '';
        robpp.calibrate = '%Values calculated by robpca.';
        robpp.out = {robpcaout.M};
        preprocessing{1} = [preprocessing{1} robpp];
        model.detail.preprocessing{1} = [model.detail.preprocessing{1} robpp];

        %copy used samples into includ field
        model.detail.includ{1} = xpp.include{1}(model.detail.robustpca.flag);
        model.detail.includ{1} = model.detail.includ{1}(:)';
        model.detail.originalinclude = {xpp.include{1}};

        %Generate SSQ table.
        sqscr = sum(model.detail.robustpca.T(model.detail.robustpca.flag,:).^2);
        model.detail.ssq(:,3) = 100*sqscr./(sum(sqscr)+sum(model.detail.robustpca.od(model.detail.robustpca.flag,:).^2));
        model.detail.ssq(:,4) = cumsum(model.detail.ssq(:,3));

        %create classes for "flagged" (excluded outliers)
        % flagged will contain 2 for user-excluded samples, 1 for robust-excluded
        % samples, and 0 for included samples.
        flagged = ones(1,size(xpp,1)).*2;  %actual excluded
        flagged(xpp.include{1}) = 1;   %ROBPCA excluded
        flagged(xpp.include{1}(model.detail.robustpca.flag)) = 0;  %included
        j = 1;
        while j<=size(model.detail.class,3);
          if isempty(model.detail.class{1,1,j});
            break;
          end
          j = j+1;
        end
        model.detail.class{1,1,j} = flagged;  %store classes there
        model.detail.classlookup{1,1,j} = {0 'Modeled'; 1 'Outlier'; 2 'User-excluded'};
        model.detail.classname{1,1,j} = 'Outlier Status';

    end
    if varargin{2}>datarank; varargin{2} = datarank; end
    
    %copy calibrated preprocessing info into model
    model.detail.preprocessing = preprocessing;   
    
    switch options.display
    case 'on'
      ssqtable(model.detail.ssq(1:varargin{2},:))
    end
        
  else
    %model from raw model (4 inputs, rawmodel=1) OR predict-from-model call (4 inputs, rawmodel=0)
    
    model   = varargin{4};
    if ~strcmp(lower(model.modeltype),'pca');
      error('Input MODEL is not a PCA model');
    end
    
    if size(x.data,2)~=model.datasource{1}.size(1,2)
      error('Variables included in data do not match variables expected by model');
    elseif length(x.include{2,1})~=length(model.detail.includ{2,1}) | any(x.includ{2,1} ~= model.detail.includ{2,1});
      missing = setdiff(model.detail.includ{2,1},x.include{2,1});
      x.data(:,missing) = nan;  %replace expected but missing data with NaN's to force replacement
      x.includ{2,1} = model.detail.includ{2,1};
    end
    
    if mdcheck(x.data(:,x.includ{2,1}));
      if strcmp(options.display,'on'); warning('Missing Data Found - Replacing with "best guess" from existing model. Results may be affected by this action.'); end
      x = replace(model,x);
    end  
    
    %apply preprocessing in model passed in
    preprocessing = model.detail.preprocessing;
    if ~isempty(preprocessing{1});
      xpp = preprocess('apply',preprocessing{1},x);
    else
      xpp = x;
    end

    %truncate loadings if # of pcs is < # of loadings
    if varargin{2} < size(model.loads{2,1},2);
      model.loads{2,1} = model.loads{2,1}(:,1:varargin{2});
    end

    if ~options.rawmodel
      model = copydsfields(x,model,1,{1 1});      %copy only mode one labels, etc.
    end

    %Update time and date.
    model.date = date;
    model.time = clock;
    
  end

  %copy options into model
  model.detail.options = options;

  if ~options.rawmodel | predictmode;    %up to here is all we need for raw model calibrate mode (predictmode=0, rawmodel=1)
    
    if options.rawmodel & predictmode; predictmode = 0; end   
    %rawmodel here means that this is really just a normal "reduce raw model" call, do NOT treat as predict
    
    %calc scores for all samples
    model.loads{1,1} = xpp.data(:,model.detail.includ{2})*model.loads{2,1};    
    
    if ~predictmode;  %use EXISTING mean/std if predict mode
%       model.detail.means{1,1}  = mean(x.data(model.detail.includ{1},:)); %mean of X-block
%       model.detail.stds{1,1}   = std(x.data(model.detail.includ{1},:));  %std of X-block
    end
    
    %calculate tsqs & residuals      
    %X-Block Statistics
    model.detail.data{1}    = x;
    model.pred{1}    = model.loads{1,1}*model.loads{2,1}';
    model.detail.res{1}     = xpp.data(:,model.detail.includ{2,1}) - model.pred{1};
    model.ssqresiduals{1,1} = model.detail.res{1}.^2;
    model.ssqresiduals{2,1} = sum(model.ssqresiduals{1,1}(model.detail.includ{1,1},:),1); %var SSQs based on cal samples only
    model.ssqresiduals{1,1} = sum(model.ssqresiduals{1,1},2);    %sample SSQs for ALL samples (excluded vars already gone)

    if ~predictmode;   %use EXISTING limit if standard predict mode
      if size(model.detail.ssq,1)>varargin{2};
        %use existing residual eigenvalues from the SSQ table if we have them (faster than calculating from raw resids)
        model.detail.reseig = model.detail.ssq(varargin{2}+1:end,2);
        if options.confidencelimit>0
          model.detail.reslim{1,1} = residuallimit(model.detail.reseig, options.confidencelimit);
        else
          model.detail.reslim{1,1} = 0;
        end
      else
        if options.confidencelimit>0
          %calculate residual eigenvalues using raw residuals matrix
          [model.detail.reslim{1,1} model.detail.reseig] = residuallimit(model.detail.res{1}(model.detail.includ{1,1},:), options.confidencelimit);
          model.detail.reseig(model.detail.reseig==0) = [];  %drop all-zero eigenvalues (not really defined)
        else
          model.detail.reslim{1,1} = 0;
        end
      end
    end
    
    if ~predictmode;   %use new scores if not predict mode
      if length(model.detail.includ{1,1})>1;
        f = sqrt(1./(diag(model.loads{1,1}(model.detail.includ{1,1},:)'* ...
          model.loads{1,1}(model.detail.includ{1,1},:))/ ...
          (length(model.detail.includ{1,1})-1)));
      else
        f = 1;
      end
    else  %use old scores and includ from original data to calculate f
      origmodel = varargin{4};
      if length(origmodel.detail.includ{1,1})>1;
        f = sqrt(1./(diag(origmodel.loads{1,1}(origmodel.detail.includ{1,1},:)'* ...
          origmodel.loads{1,1}(origmodel.detail.includ{1,1},:))/ ...
          (length(origmodel.detail.includ{1,1})-1)));
      else
        f = 1;
      end
    end
    model.tsqs{1,1}          = sum((model.loads{1,1}*diag(f)).^2,2);
    model.tsqs{2,1}          = sum(model.loads{2,1}.^2,2)'*(length(model.detail.includ{2,1})-1);

    if ~predictmode;   %use EXISTING limit if standard predict mode
      if length(model.detail.includ{1,1})>varargin{2} & options.confidencelimit>0
        model.detail.tsqlim{1,1} = tsqlim(length(model.detail.includ{1,1}),varargin{2},options.confidencelimit*100);
      else
        model.detail.tsqlim{1,1} = 0;
      end    
    end
        
    if predictmode
      model.detail.rmsec  = [];
      model.detail.rmsecv = [];
      model.detail.cvpred = [];
      model.detail.cv     = '';
      model.detail.split  = [];
      model.detail.iter   = [];
    end
    
    %test for Compact Blockdetails with calibrate mode
    if ~predictmode & strcmp(options.blockdetails,'compact')
      warning(['Compact model not yet available - will return standard model'])
      options.blockdetails = 'standard';
    end
    
    %handle model compression
    switch lower(options.blockdetails)
      case 'compact'
        if predictmode
          pred = modelstruct(model.modeltype,1);  %get reduced model structure
          pred.datasource = model.datasource;
          pred.date       = model.date;
          pred.time       = model.time;
          pred.scores     = model.loads{1};
          pred.pred       = model.pred{1};
          pred.tsqs       = model.tsqs{1};
          pred.ssqresiduals = model.ssqresiduals{1};
          model = pred;
        else
          %NOTE - developmental compact model code - test above locks this out
          model.info = ['Compact ' model.description{1}];
          %keep these fields, remove all others
          for field = setdiff(fieldnames(model)',{ 'modeltype' 'datasource' 'date' 'time' 'info' 'loads' 'detail' });
            model = rmfield(model,field{1});
          end
          for field = { 'includ' 'preprocessing' };
            model = setfield(model,field{1},getfield(model.detail,field{1}));
          end
        end
      case 'standard'
        if predictmode
          model.modeltype = [model.modeltype '_PRED'];
        end
        model.detail.data{1} = [];
        model.pred{1} = [];
        model.detail.res{1}  = [];
    end
    
  end
  
  %output version options  
  switch options.outputversion
  case 2
    %I/O: [scores,loads,ssq,res,reslm,tsqlm,tsq] = pca(data,plots,scl,lvs);
    varargout{1} = model.loads{1}; 
    varargout{2} = model.loads{2}; 
    varargout{3} = model.detail.ssq; 
    varargout{4} = model.ssqresiduals{1}; 
    varargout{5} = model.detail.reslim{1}; 
    varargout{6} = model.detail.tsqlim{1};
    varargout{7} = model.tsqs{1};
  otherwise
    varargout{1} = model;
  end  

  try
    switch lower(options.plots)
      case 'final'
        if ~predictmode
          plotloads(model);
          plotscores(model);
        else
          plotloads(varargin{4},model);
          plotscores(varargin{4},model);
        end
    end
  catch
    warning(lasterr)
  end
end
%End Input

%--------------------------
function out = optiondefs()

defs = {
  
%name                    tab              datatype        valid                            userlevel       description
'display'                'Display'        'select'        {'on' 'off'}                     'novice'        'Governs level of display.';
'plots'                  'Display'        'select'        {'none' 'final'}                 'novice'        'Governs level of plotting.';
'outputversion'          'Standard'       'select'        {2 3}                            'novice'        'Governs output format.';
'algorithm'              'Standard'       'select'        {'svd' 'maf' 'robustpca'}        'novice'        'Algorithm for decomposition. ''svd'' is standard decomposition algorithm. ''robustpca'' uses the robust PCA algorithm of the LIBRA toolbox (automatic outlier exclusion). ''maf'' is Maximum Autocorrelative Factors and requires Eigenvector''s MIA_Toolbox.';
'preprocessing'          'Standard'       'matrix'        ''                               'novice'        'Preprocessing structure.';
'blockdetails'           'Standard'       'select'        {'compact' 'standard' 'all'}     'novice'        'Extent of predictions and raw residuals included in model. ''standard'' = none, ''all'' x-block.';
'confidencelimit'        'Standard'       'double'        'float(0:1)'                     'novice'        'Confidence level for Q and T2 limits (fraction between 0 and 1)';
'roptions.alpha'         'Robust Options' 'double'        'float(.5:1)'                     'novice'        '(1-alpha) measures the number of outliers the algorithm should resist. Any value between 0.5 and 1 may be specified (default = 0.75). These options are only used when algorithm is robust PCA.';
};

out = makesubops(defs);

