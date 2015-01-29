function [ax,mx,stdx,msg] = auto(x,options)
%AUTO Autoscales matrix to mean zero unit variance.
%  INPUT:
%        x = MxN matrix to autoscale.
%
%  OPTIONAL INPUT:
%   options = structure array with the following fields:
%              offset: scales by stdx+offset {default = 0}:
%                      if (offset = inf) no scaling used (same as MNCN),
%                      if (offset = scalar), this scalar offset is
%                      applied to all N columns, or
%                      (offset = 1xN vector) allows column specific
%                      scaling.
%             display: [ {'off'} | 'on'] governs level of display.
%    matrix_threshold: scalar threshold of fraction of missing data
%                      in (x) {default = 0.15}, and
%    column_threshold: scalar threshold of fraction of missing data
%                      in a single column {default = 0.25}.
%           algorithm: [ {'standard'} | 'robust'] scaling algorithm.
%                      'robust' uses MADC for scaling and median instead of
%                      mean. Should be used for robust techniques.
%        stdthreshold: [0] scalar or vector of standard deviation threshold
%                      values. If a standard deviation is below its
%                      corresponding threshold value, the threshold value
%                      will be used in lieu of the actual value. Note that
%                      the actual standard deviation is always returned,
%                      whether or not it exceedes the threshold. A scalar
%                      value is used as a threshold for all variables.
%      badreplacement: [0] value to use in place of standard deviation
%                      values of 0 (zero). Typical values used with the
%                      following effects:
%                        0 = Any value in given variable is set to zero.
%                            Variable is effectively excluded (but still
%                            expected by model). This is also the behavior
%                            when badreplacement = inf.
%                        1 = Values different from mean of the given
%                            variable are flagged in Q residuals with no
%                            reweighting.
%                        Values >0 and <inf give the variable different
%                            weighting in the Q residuals (values >1
%                            down-weight the bad variables for Q residual
%                            calculations, values <1 up-weight the bad
%                            variables.)
%
%  OUTPUTS:
%       ax = MxN scaled and centered matrix, e.g.,
%            with mean-zero unit variance columns.
%       mx = 1xN vector of means (mx) used for centering.
%     stdx = 1xN vector of standard deviations (stdx) used for scaling.
%      msg = returns any warning messages.
%  If missing data (NaNs) are found, the available data are autoscaled if
%  the fraction missing is not above the thresholds specified in (options).
%
%I/O: [ax,mx,stdx,msg] = auto(x,options);
%I/O: [ax,mx,stdx,msg] = auto(x,offset);
%I/O: auto demo
%
%See also: GSCALE, MEDCN, MNCN, NORMALIZ, NPREPROCESS, REGCON, RESCALE, SCALE, SNV

% Copyright © Eigenvector Research, Inc. 1991-2008
% Licensee shall not re-compile, translate or convert "M-files" contained
%  in PLS_Toolbox for use with any software other than MATLAB®, without
%  written permission from Eigenvector Research, Inc.
%Modified 11/93
%Checked on MATLAB 5 by BMW  1/4/97
%Added support for non-doubles (do by cols) by JMS 6/21/2001
%Added offset support JMS 9/18/01
%JMS 6/6/02 -Incorporated missing data support
%  -standard IO mods
%JMS 2/7/03 -Do missing data calculation on ONLY those columns that need it
%JMS 9/29/03 -Modified logic to detect zero-variance columns and method of
%  zeroing them out
%JMS 10/27/03 -Fixed zero-norm bug (introduced by previous fix)
%JMS 1/23/04 -Added message warning if zero-standard deviation column is found
%  -copied most robust standard deviation code throughout file
%JMS 5/10/04 -fixed one-row bug

if nargin == 0
    x = 'io'; 
end

varargin{1} = x;
if ischar(varargin{1});
    
  options = [];
  options.name = 'options';
  options.offset = 0;
  options.display = 'off';
  options.matrix_threshold = .15;    %Thresholds for "too much missing data"
  options.column_threshold = .25;
  options.algorithm = 'standard'; %Algorithm for determining
  options.stdthreshold = 0;     %no threshold
  options.badreplacement = 0;   %replace bad stds with this value

%   if nargout==0; 
%       evriio(mfilename,varargin{1},options); 
%   else 
%       ax = evriio(mfilename,varargin{1},options); 
%   end
  
  ax = options;
  
  return;
  
end

if nargin < 2
  options = auto('options');
else
  if isa(options,'double');
    offset         = options;
    options        = [];
    options.offset = offset;
  end
  options = reconopts(options,'auto');
end

msg     = {};

%check for DataSet object
originaldso = [];
if isdataset(x)
  originaldso = x;
  incl = x.include;
  if length(options.stdthreshold)==size(x,2);
    %apply include field to stdthreshold (if applicable)
    options.stdthreshold = options.stdthreshold(incl{2});
  end
  %extract data from x for initial calculation.
  x = originaldso.data(incl{:});
end

%check for missing data
pctmd = mdcheck(x);
if pctmd > options.matrix_threshold;
  error(['Too much missing data to analyze'])
end
if pctmd > 0
  msg{end+1} = 'Missing data found. Using missing data logic.';
  if strcmp(options.display,'on')
    disp(msg{end});
  end
end

[m,n] = size(x);

%Check the threshold for validity
if isempty(options.stdthreshold)
  options.stdthreshold = 0;
end
%test length of threshold vector compared to data
if length(options.stdthreshold)>1 && length(options.stdthreshold)~=n
  error('Standard Deviation Threshold does not match number of included variables')
end
%make sure it is a vector (auto-expands scalar values to vector)
options.stdthreshold = ones(1,n).*options.stdthreshold;

wrn=warning;
warning off;
if isa(x,'double')

  %first pass - do all columns (will be recalc'ed for columns w/missing data)
  switch options.algorithm
    case 'standard'
      mx    = mean(x,1);
    case 'robust'
      mx    = median(x,1);
  end
  if isfinite(options.offset)
    switch options.algorithm
      case 'standard'
        stdx          = std(x,[],1)+options.offset;
      case 'robust'
        stdx          = madc(x)+options.offset;
    end

    %test against threshold
    stdxapply  = max([stdx;options.stdthreshold],[],1);
    
    %look for bad standard deviations
    badstdx    = findzerostd(stdxapply,x);
    if options.badreplacement==0;
      stdx(badstdx)      = 0;
      stdxapply(badstdx) = inf;
    else
      stdxapply(badstdx) = options.badreplacement;
      stdx(badstdx)      = options.badreplacement;
    end
    
    if any(badstdx);
      for loop = find(badstdx);
        msg{end+1} = ['Column ' num2str(loop) ' has a standard deviation of zero.'];
        if strcmp(options.display,'on')
          disp(msg{end});
        end
      end
    end
  else
    %infinite offset, no scaling at all
    stdx          = ones(1,n);
    stdxapply     = stdx;
  end
  rvect = ones(m,1);
  ax    = (x-rvect*mx)./(rvect*stdxapply);
  
  if pctmd > 0;
    %redo columns with missing data
    [junk,mxmdmap] = mdcheck(mx);
    for loop = find(mxmdmap);
      onecolumn    = x(:,loop);
      [pctmd_c,mapmd_c] = mdcheck(onecolumn);
      if pctmd_c == 0; mapmd_c = zeros(m,1); end;
      if pctmd_c > options.column_threshold;
        error(['Too much missing data in column ' num2str(loop) ' to analyze'])
      end
      switch options.algorithm
        case 'standard'
          mx(1,loop)   = mean(onecolumn(~mapmd_c));
        case 'robust'
          mx(1,loop)   = median(onecolumn(~mapmd_c));
      end
      if isfinite(options.offset)
        %stdx(1,loop)  = std(onecolumn(~mapmd_c))+options.offset;
        switch options.algorithm
          case 'standard'
            stdx(1,loop)  = std(onecolumn(~mapmd_c))+options.offset;
          case 'robust'
            stdx(1,loop)  = madc(onecolumn(~mapmd_c))+options.offset;
        end

        %test against threshold
        stdxapply(1,loop) = max(stdx(1,loop),options.stdthreshold(1,loop));

        %check for "bad" (zero) values
        if findzerostd(stdxapply(1,loop),onecolumn(~mapmd_c));
          if options.badreplacement==0;
            stdx(1,loop)      = 0;
            stdxapply(1,loop) = inf;
          else
            stdx(1,loop)      = options.badreplacement;
            stdxapply(1,loop) = options.badreplacement;
          end
          msg{end+1} = ['Column ' num2str(loop) ' has a standard deviation of zero.'];
          if strcmp(options.display,'on')
            disp(msg{end});
          end
        end
      else
        %infinite offset, no scaling at all
        stdx(1,loop) = 1;
        stdxapply(1,loop) = stdx(1,loop);
      end
      ax(:,loop)   = (onecolumn-mx(1,loop))./stdxapply(1,loop);
    end
  end

else  %other than double, convert one column at a time

  ax = feval(class(x),zeros(m,n));
  for loop = 1:n;
    onecolumn    = double(x(:,loop));
    [pctmd_c,mapmd_c] = mdcheck(onecolumn);
    if pctmd_c == 0; mapmd_c = zeros(m,1); end;
    if pctmd_c > options.column_threshold;
      error(['Too much missing data in column ' num2str(loop) ' to analyze'])
    end
    switch options.algorithm
      case 'standard'
        mx(1,loop)   = mean(onecolumn(~mapmd_c));
      case 'robust'
        mx(1,loop)   = median(onecolumn(~mapmd_c));
    end
    
    if isfinite(options.offset)
      %stdx(1,loop)  = std(onecolumn(~mapmd_c))+options.offset;

      switch options.algorithm
        case 'standard'
          stdx(1,loop)  = std(onecolumn(~mapmd_c))+options.offset;
        case 'robust'
          stdx(1,loop)  = madc(onecolumn(~mapmd_c))+options.offset;
      end

      %test against threshold
      stdxapply(1,loop) = max(stdx(1,loop),options.stdthreshold(1,loop));

      %check for "bad" (zero) values
      if findzerostd(stdxapply(1,loop),onecolumn(~mapmd_c));
        if options.badreplacement==0;
          stdx(1,loop)      = 0;
          stdxapply(1,loop) = inf;
        else
          stdx(1,loop)      = options.badreplacement;
          stdxapply(1,loop) = options.badreplacement;
        end
        msg{end+1} = ['Column ' num2str(loop) ' has a standard deviation of zero.'];
        if strcmp(options.display,'on')
          disp(msg{end});
        end
      end

    else
      stdx(1,loop) = 1;
      stdxapply(1,loop) = stdx(1,loop);
    end
    ax(:,loop)   = feval(class(x),(onecolumn-mx(1,loop))./stdxapply(1,loop));
  end

end

if isdataset(originaldso);
  %if we started with a DSO, re-insert back into DSO
  scaleddata = originaldso.data*nan;  %block out all data
  scaleddata(incl{:}) = ax;  %insert scaled data (included columns only)
  originaldso.data = scaleddata;
  ax = originaldso;
end

msg = char(msg);

warning(wrn);

%--------------------------------------------
function badstdx = findzerostd(stdx,x)

nrm           = sqrt(sum(x.^2,1));
nrm(nrm==0)   = 1/eps;
badstdx       = (stdx./nrm)<eps*3;
