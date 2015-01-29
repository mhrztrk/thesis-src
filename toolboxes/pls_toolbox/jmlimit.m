function rescl = jmlimit(pc,s,cl)
%JMLIMIT Confidence limits for Q residuals via Jackson-Mudholkar.
%  Inputs are the number of PCs used (pc), the vector of 
%  eigenvalues (s), and the confidence limit (cl) expressed
%  as a fraction (e.g. 0.95). The output (rescl) is the confidence limit
%  based on the method of Jackson and Mudholkar. See CHILIMIT for an 
%  alternate method of residual limit calculation based on chi^2.
%
%I/O: rescl = jmlimit(pc,s,cl);
%
%Example: rescl = jmlimit(2,ssq(:,2),0.95);
%         rescl = jmlimit(4,model.detail.ssq(:,2),0.99); %for PCA a model
%
%See also: ANALYSIS, CHILIMIT, PCA, RESIDUALLIMIT

%Licensee shall not re-compile, translate or convert "M-files" contained
% in PLS_Toolbox for use with any software other than MATLAB®, without
% written permission from Eigenvector Research, Inc.
%Copyright Eigenvector Research, Inc. 1997-2008
%nbg 4/97,7/97
%bmw 12/99
%jms 4/03 -revised to allow multiple limits to be calculated in one call
%  (vector of cl)

if nargin == 0; pc = 'io'; end
varargin{1} = pc;
if ischar(varargin{1});
  options = [];
  if nargout==0; clear rescl; evriio(mfilename,varargin{1},options); else; rescl = evriio(mfilename,varargin{1},options); end
  return; 
end

if cl>=1|cl<=0
  error('confidence limit must be 0<cl<1')
end
cl = cl*100;
[m,n] = size(s);
if m>1&n>1
  error('input s must be a vector')
end
if n>m
  s   = s';
  m   = n;
end
if pc>=length(s)
  rescl = 0;
else
  cl     = 2*cl-100;
  theta1 = sum(s(pc+1:m,1));
  theta2 = sum(s(pc+1:m,1).^2);
  theta3 = sum(s(pc+1:m,1).^3);
  if theta1==0;
    rescl = 0;
  else
    h0     = 1-2*theta1*theta3/3/(theta2.^2);
    if h0<0.001
      h0 = 0.001;
    end
    ca    = sqrt(2)*erfinv(cl/100);
    h1    = ca*sqrt(2*theta2*h0.^2)/theta1;
    h2    = theta2*h0*(h0-1)/(theta1.^2);
    rescl = theta1*(1+h1+h2).^(1/h0);
  end
end
