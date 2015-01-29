function fstat = ftest(p,n,d,flag)
%FTEST F test and inverse F test statistic.
%  For (flag) set to 1 {default} FTEST calculates the
%  F statistic (fstat) given the probability point (p)
%  and the numerator (n) and denominator (d) degrees
%  of freedom. (flag is an optional input {default = 1}.
%
%  For (flag) set to 2 FTEST calculates the probability
%  point (fstat) given the F statistic (p) and the
%  numerator (n) and denominator (d) degrees of freedom.
%
%Example:
%  a = ftest(0.05,5,8)
%    a = 3.6875
%  a = ftest(3.6875,5,8,2)
%    a = 0.050
%  
%I/O: fstat = ftest(p,n,d,flag);
%I/O: ftest demo
%
%See also: ANOVA1W, ANOVA2W, STATDEMO, TTESTP

%Copyright Eigenvector Research, Inc. 1991-2008
%Licensee shall not re-compile, translate or convert "M-files" contained
% in PLS_Toolbox for use with any software other than MATLAB®, without
% written permission from Eigenvector Research, Inc.
%the orginal was based on a public domain stats toolbox
%Modified 11/93,12/94 BMW, 10/96,12/97 NBG
%modified nbg 3/02: error trapping, switch/case, added tol, change to golden section
% numbers are better behaved, changed help

if nargin == 0; p = 'io'; end
if ischar(p);
  options = [];
  if nargout==0; evriio(mfilename,p,options); else; fstat = evriio(mfilename,p,options); end
  return; 
end

if prod(size(n))>1 | prod(size(d))>1           %nbg 3/02
  error('Inputs N and D must be scalar.')
end
n    = n/2;
d    = d/2;
if nargin<4, flag = 1; end

switch flag
case 1
  p    = 1-p;                 %e.g. change 0.05 to 0.95
  ic   = 1;     itmax = 30;   %iteration count, max it count 3/02
  xtol = 1e-8;  ftol = 1e-6;  %tolerance on x,  tol on function betainc(x)-p  3/02
  xl   = 0.0;                 %x left
  xr   = 1.0;                 %x right
  fxl  = -p;                  %betainc(0)-p
  fxr  = 1.0 - p;             %betainc(1)-p
  if fxl*fxr>0
    error('probability not in the range(0,1) ')
  else
%     while ic<itmax
%       x   = (xl+xr)*0.5;
%       p1  = betainc(x,n,d);
%       fcs = p1-p;
%       if fcs*fxl>0
%         xl  = x;
%         fxl = fcs;
%       else
%         xr  = x;
%         fxr = fcs;
%       end
%       xrmxl = xr-xl;
%       if xrmxl<=xtol | abs(fcs)<=ftol %3/02 changed to xtol,ftol reduced ftol from 1e-4 to 1e-6
%         break
%       else
%         ic  = ic+1;
%       end
%     end
    a     = 0.382;
    while ic<itmax   %based on golden section search, might use Newton-Raphson
      x   = xl+(xr-xl)*a;    %(xl+xr)*0.5;
      fcs = betainc(x,n,d)-p;
      if fcs*fxl>0
        xl  = x;
        fxl = fcs;
      else
        xr  = x;
        fxr = fcs;
      end
      xrmxl = xr-xl;
      if any([xrmxl<=xtol abs(fcs)<=ftol]) %3/02 changed to xtol,ftol reduced ftol from 1e-4 to 1e-6
        break
      else
        ic  = ic+1;
      end
    end
  end
  if ic>=itmax
    warning(['FTEST failed to converge. Number of iterations at maximum ',int2str(itmax)])
  end
  % Inverted numerator and denominator 12-26-94
  fstat = (d * x) / (n - n * x);
case 2
  fstat = betainc(2*d./(2*d+2*n*p),d,n);  %03/02 changed / to ./ to allow for vectors
otherwise
  error('Input FLAG not recognized.')
end
