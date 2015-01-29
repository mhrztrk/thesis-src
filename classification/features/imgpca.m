function model = imgpca(mwa,scaling,nocomp)
%IMGPCA Principal Components Analysis of Multivariate Images.
%  IMGPCA uses principal components analysis to make
%  psuedocolor maps of multivariate images. The input is the
%  multivariate image (mwa). Optional inputs are (scaling) the
%  scaling to be used, and the number of PCs to calculate (nocomp).
%    scaling = 'auto' uses autoscaling {default},
%    scaling = 'mncn' uses mean centering, and
%    scaling = 'none' uses no scaling.
%
%  It is assumed that the image (mwa) is a 3 dimensional (m x n x p)
%  array where each image is m x n pixels and there are p images.
%  IMGPCA presents each scores, residual, and T^2 matrix as a
%  psuedocolor image. If 3 are more PCs are selected (nocomp>=3),
%  a composite of the first three  PCs is shown as an rgb image,
%  with red for the first PC, green for the second, and blue for the
%  the third.
%
%  The output (model) is a structure with the following fields:
%
%     xname: input data name
%      name: type of model, always 'IPCA'
%      date: date of model creation
%      time: time of model creation
%      size: dimensions of input data
%    nocomp: number of PCs in model
%     scale: type of scaling used
%     means: mean vector for PCA model
%      stds: standard deviation vector for PCA model
%       ssq: variance captured table data
%    scores: PCA scores stored as m x n x nocomp array (uint8)
%     range: original range of PCA scores before mapping to uint8
%     loads: PCA loadings
%       res: PCA residuals stored as m x n array (uint8)
%    reslim: Q limit
%       tsq: PCA T^2 values stared as m x n array (unit8)
%    tsqlim: T^2 limit
%
%  Note that the scores, residuals and T^2 matrices are stored
%  as unsigned 8 bit integers (uint8) scaled so their range is 
%  0 to 255. These can be viewed with the IMAGE function, but 
%  be sure the current colormap has 256 colors. For example, to
%  view the scores on the second PC using the jet colormap:
%
%   image(model.scores(:,:,2)), colormap(jet(256)), colorbar
%
%I/O: model = imgpca(mwa,scaling,nocomp);
% 
% IMGPCA can also be used to apply existing IPCA models to 
% new images as follows:
%
%I/O: newmod = imgpca(mwa,model,plots);
%
% If plots == 0, no plots are produced.
%
%See also: CONTRASTMOD, IMAGEGUI, IMGSELCT, IMGSIMCA, IMREAD, ISIMCAPR

%Copyright Eigenvector Research, Inc. 1998-2004
%Licensee shall not re-compile, translate or convert "M-files" contained
% in PLS_Toolbox for use with any software other than MATLAB®, without
% written permission from Eigenvector Research, Inc.
%BMW
%nbg 11/00 changed help
%BMW 9/02 made to accept DSOs, evriio
%BNW 12/02 added labels to loads plots
%
% code modified by K.Artyushkova not to display the % variance captured and not to display
% figures

    if nargin == 0; mwa = 'io'; end
    varargin{1} = mwa;
    if ischar(varargin{1});
      options = [];
      if nargout==0 
          evriio(mfilename,varargin{1},options); 
      else
          model = evriio(mfilename,varargin{1},options); 
      end
      return; 
    end

    smwa = size(mwa);
    if isa(mwa,'dataset')
      inds = mwa.includ;
      s = ' ';
      varlabel = [s(ones(length(inds{3}),1)) mwa.label{3}(inds{3},:)];
      mwa = mwa.data(inds{:});
    else
      varlabel = [];
    end

    % If scaling not specified, set it to auto scaling
    if (nargin < 2 || isempty(scaling))
      scaling = 'auto';
    end

    if isa(scaling,'char') % Not a model, so make one.

    ms = size(mwa);
    nr = ms(1)*ms(2);
    mwa = reshape(mwa,nr,ms(3));

    if strcmp(class(mwa),'uint8') || strcmp(class(mwa),'double')

        % Calculate the scatter matrix
        scmat = zeros(ms(3),ms(3));
        for i = 1:ms(3)
            for j = 1:i  
                scmat(i,j) = double(mwa(:,i))'*double(mwa(:,j));
                if i ~= j
                    scmat(j,i) = scmat(i,j);
                end
            end
        end

        % Scale data as specified
        if (strcmp(scaling,'auto'))
            scaling = 'auto';
            mns = mean(mwa);
            stds = sqrt((diag(scmat)' + nr*mns.^2 - 2*sum(mwa).*mns)/(nr-1)); 
            scmat = (inv(diag(stds))*(scmat - mns'*mns*nr)*inv(diag(stds)))/(nr-1);
        elseif(strcmp(scaling,'mncn'))
            scaling = 'mncn';
            mns = mean(mwa);
            scmat = (scmat - mns'*mns*nr)/(nr-1);
            stds = ones(1,ms(3));
        elseif (strcmp(scaling,'none'))
            scaling = 'none';
            scmat = scmat/(nr-1);
            mns = zeros(1,ms(3));
            stds = ones(1,ms(3));
        else
            error('Scaling not of known type')
        end

        % Calculate the loadings
        if nargin < 3
            [~,s,loads] = svd(scmat);
            nocomp = ms(3);
        else
            [~,s,loads] = svd(scmat);
        end

        % Change the sign on the loads to be mostly positive
        loads = loads*diag(sign(sum(loads)));

        % Display the variance captured table.
        if strcmp(scaling,'none')
            disp('  ')
            disp('Warning: Data was not mean centered.')
            disp(' Variance captured table should be read as sum of')
            disp(' squares captured.') 
        end

        temp = diag(s(1:nocomp,1:nocomp))*100/(sum(diag(s)));
        ssq  = [(1:nocomp)' diag(s(1:nocomp,1:nocomp)) temp cumsum(temp)];
        disp('   ')
        disp('        Percent Variance Captured by PCA Model')
        disp('  ')
        disp('Principal     Eigenvalue     % Variance     % Variance')
        disp('Component         of          Captured       Captured')
        disp(' Number         Cov(X)        This  PC        Total')
        disp('---------     ----------     ----------     ----------')
        format = '   %3.0f         %3.2e        %6.2f         %6.2f\n';
        mprint = min([20 nocomp]);
        for i = 1:mprint
            fprintf(1,format,ssq(i,:));
        end
        
        fprintf(1,'\n');
        for i=1:nocomp
            fprintf(1,'PC%d = ', i);
            for j=1:ms(3)
                fprintf(1,'%7.3f * mwa(%d)', loads(j,i), j);
                if(j ~= ms(3))
                    fprintf(1,' + ');
                end
            end
            fprintf(1,'\n');
        end
        fprintf(1,'\n');
        
        % Calculate the scores
        scores = uint8(zeros(nr,nocomp));
        scr = zeros(2,nocomp+2);
        for j = 1:nocomp
            ts = zeros(nr,1);
            if strcmp(scaling,'none')
%                 for i = 1:nr
%                     ts(i,:) = double(mwa(i,:))*loads(:,j);
%                 end
                A = bsxfun(@times, double(A), loads(:,j)');
                ts(i,:) = sum(A,2);
                
            elseif strcmp(scaling,'mncn')            
%                 for i = 1:nr
%                     ts(i,:) = (double(mwa(i,:))-mns)*loads(:,j);
%                 end
                
                A = bsxfun(@minus, double(mwa), mns);
                A = bsxfun(@times, double(A), loads(:,j)');
                ts(i,:) = sum(A,2);
                
            elseif strcmp(scaling,'auto')           
%                 for i = 1:nr
%                    ts(i,:) = ((double(mwa(i,:))-mns)./stds)*loads(:,j);
%                 end  
                A = bsxfun(@minus, double(mwa), mns);
                A = bsxfun(@rdivide, A, stds);
                A = bsxfun(@times, double(A), loads(:,j)');
                ts = sum(A,2);
                
            end
            scr(1,j) = min(ts);
            scr(2,j) = max(ts);
            scores(:,j) = round(255*(ts-min(ts))/max(ts-min(ts)));
        end
        clear A ts
        
        % Calculate the residuals and T^2
        imppt = eye(ms(3))-loads*loads';
        res = zeros(nr,1);
        tsqs = zeros(nr,1);
        if strcmp(scaling,'none')
            for i = 1:nr
                smwa = double(mwa(i,:));
                res(i) = sum((smwa*imppt).^2);
                tsqs(i) = sum(((smwa*loads)./sqrt(ssq(1:nocomp,2)')).^2);
            end
        elseif strcmp(scaling,'mncn')
            for i = 1:nr
                smwa = double(mwa(i,:))-mns;
                res(i) = sum((smwa*imppt).^2);
                tsqs(i) = sum(((smwa*loads)./sqrt(ssq(1:nocomp,2)')).^2);
            end
        elseif strcmp(scaling,'auto')
            for i = 1:nr
                smwa = (double(mwa(i,:))-mns)./stds;
                res(i) = sum((smwa*imppt).^2);
                tsqs(i) = sum(((smwa*loads)./sqrt(ssq(1:nocomp,2)')).^2);
            end
        end
        scr(1,nocomp+1) = min(res);
        scr(2,nocomp+1) = max(res);
        res = uint8(255*res/max(res));
        scr(1,nocomp+2) = min(tsqs);
        scr(2,nocomp+2) = max(tsqs);
        tsqs = uint8(255*tsqs/max(tsqs));

        % Calculate the Q limit
        if nocomp < ms(3);
            temp = diag(s);
            emod = temp(nocomp+1:end);
            th1 = sum(emod);
            th2 = sum(emod.^2);
            th3 = sum(emod.^3);
            h0 = 1 - ((2*th1*th3)/(3*th2^2));
            if h0 <= 0.0
                h0 = .0001;
                disp('  ')
                disp('Warning:  Distribution of unused eigenvalues indicates that')
                disp('          you should probably retain more PCs in the model.')
            end
            q = th1*(((1.65*sqrt(2*th2*h0^2)/th1) + 1 + th2*h0*(h0-1)/th1^2)^(1/h0));
            disp('  ')
            str = sprintf('The 95 Percent Q limit is %g',q);
            disp(str)
        else
            q = 0;
        end

        %  Calculate T^2 limit using ftest routine
        if nr > 300
            tsq = (nocomp*(nr-1)/(nr-nocomp))*ftest(.05,nocomp,300);
        else
            tsq = (nocomp*(nr-1)/(nr-nocomp))*ftest(.05,nocomp,nr-nocomp);
        end
        disp('  '), str = sprintf('The 95 Percent T^2 limit is %g',tsq); disp(str)
    end

    % Fold the scores, residuals and T^2s back up
    scores = reshape(scores,ms(1),ms(2),nocomp);
    res = reshape(res,ms(1),ms(2));
    tsqs = reshape(tsqs,ms(1),ms(2));

    model = struct('xname',inputname(1),'name','IPCA','date',date,'time',clock,...
      'size',ms,'nocomp',nocomp,'scale',scaling,'means',mns,'stds',stds,...
      'ssq',ssq,'scores',scores,'range',scr,'loads',loads,'res',res,'reslim',q,...
      'tsq',tsqs,'tsqlim',tsq);

    end
end
