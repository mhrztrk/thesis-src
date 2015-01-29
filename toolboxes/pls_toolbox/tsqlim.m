function tsqcl = tsqlim(m,pc,cl)
%TSQLIM Confidence limits for Hotelling's T^2.
%  Inputs are the number of samples (m), the number of PCs used (pc), and
%  the confidence limit (cl) where 0 < cl < 1. Optionally, (m) and (pc) can
%  be omitted and a standard model structure can be passed (model) along
%  with just the desired confidence limit (cl).
%
%  The output (tsqcl) is the confidence limit.
%
%Examples: tsqcl = tsqlim(15,2,0.95);
%          tsqcl = tsqlim(mymodel,0.95);
%
%I/O: tsqcl = tsqlim(m,pc,cl);
%I/O: tsqcl = tsqlim(model,cl);
%I/O: tsqlim demo
%
%See also: ANALYSIS, PCA, PCR, PLS

%Copyright Eigenvector Research, Inc. 1997-2008
%Licensee shall not re-compile, translate or convert "M-files" contained
% in PLS_Toolbox for use with any software other than MATLAB®, without
% written permission from Eigenvector Research, Inc.
%nbg 4/97
%nbg 3/22/02 additional error checking
%jms 6/23/03 allow multiple limits for simultaneous detn.
%jms 11/29/05 allow passing of standard model structure

if nargin == 0; m = 'io'; end
varargin{1} = m;
if ischar(varargin{1});
  options = [];
  if nargout==0; clear tsqcl; evriio(mfilename,varargin{1},options); else; tsqcl = evriio(mfilename,varargin{1},options); end
  return; 
end

if isfield(m,'modeltype')
  if nargin<2
    error('TSQLIM requires a confidence limit along with the model');
  end
  %(model,cl)
  cl = pc;
  pc = size(m.loads{2,1},2);  %number of components
  m  = length(m.detail.includ{1});
elseif nargin<3
  error('TSQLIM requires 3 inputs.')
end

tsqcl_all = [];
for cl = cl;
  
  %assume that, typical CL will be ~0.9 to 0.999
  if cl<=0
    error('INPUT CO must be >0 (0<cl<1).')
  elseif cl>=100
    error('INPUT CO must be <1 (0<cl<1).')
  elseif cl>0&cl<1
    cl = cl*100;  %convert to %
  end
  
  warning off backtrace
  if m==pc
    tsqcl = NaN;
    warning('Warning: TSQLIM not defined when pc==m. tsqcl = NaN')
  elseif pc>m
    tsqcl = NaN;
    warning('Warning: TSQLIM not defined when pc>m. tsqcl = NaN')
  else
    alpha = (100-cl)/100;
    tsqcl = pc*(m-1)/(m-pc)*ftest(alpha,pc,m-pc);
  end
  tsqcl_all(end+1) = tsqcl;
end

tsqcl = tsqcl_all;

