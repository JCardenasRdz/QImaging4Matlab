function [MAPS,MAPScfit]=qMAPS(Images_Matrix,TR,mask)
% Process imaging data with an arbiritary function
%
% [MAPS,MAPScfit]=MAPSMap(Images_Matrix,TR,mask)
%%  Inputs
%   Images_Matrix=    3D Matrix of MR Images of size (N,M,K)
%   TR=               Vector of values for xdata in LSQcurvefit of size K
%   mask=             Binariy Mask to process only an specific ROI
%
%%  OutPuts
%   MAPS.Map=  Map for the T1 of size (N,M);
%   MAPS.rsq=  Map for the Goodness of fit of size (N,M);
%
%%  Author
%   Julio Cárdenas-Rodríguez
%       The University of Arizona
%       cardenaj@email.arizona.edu
%%  Algorithm
%1) Define options for lsqcurvefit method
options = optimoptions('lsqcurvefit');
options.MaxIter=1E3;
options.MaxFunEvals=2E3;
options.TolX=1e-6;             options.TolFun=1e-6;
options.Display='off';

MAPSfunc=@(pars,xdata) pars(1).*(1-exp(-xdata./pars(2)) );

%  2) Reshape mask to vectorize and preallocate maps
[rows,cols,~]=size(Images_Matrix);
mask=reshape(mask,rows*cols,[]) ;
MAPS.map=nan(size(mask));
MAPS.rsq=MAPS.map;

%  3)  Find indexes for pixeles to be analyzed
indices= find(mask);

%  4)  preallocate MAPScfit cell
MAPScfit=cell(size(indices));

%  5)  Reshape T1 vTR 3D matrix into a 2D matrix
Images_Matrix=reshape(Images_Matrix,rows*cols,[]) ;

%  6) CurveFitting for each voxel
x0=[1.2,3.0];
lb=[1,1];
ub=[2,5];
N=length(TR);

for q=1:length(indices);
    
    Signal= Images_Matrix (indices(q),:)';   Signal=Signal./max(Signal);

[beta,~,resid,~,~,~,J]=lsqcurvefit(MAPSfunc,...
                                              x0,TR,Signal,lb,ub,options);
    
    T1    = beta(2);
    T1ci  = nlparci(beta,resid,'jacobian',J);
    rmse  =rsquare(Signal,  MAPSfunc(beta,TR)   );
    
    MAPScfit{q}.cfit.T1=     T1;
    MAPScfit{q}.cfit.T1ci=   T1ci;
    MAPScfit{q}.cfit.MAPS.map=   MAPS.map;
    
    MAPS.rsq(indices(q),1)=rmse;
    MAPS.map(indices(q),1)=T1;
end

%  Reshapes MAPS back to square matrices

MAPS.map=reshape(MAPS.map,rows,cols); 
MAPS.rsq=reshape(MAPS.rsq,rows,cols);

%%%
function [r2, rmse] = rsquare(y,f,varargin)
% Compute coefficient of determination of data fit model and RMSE
%
% [r2 rmse] = rsquare(y,f)
% [r2 rmse] = rsquare(y,f,c)
%
% RSQUARE computes the coefficient of determination (R-square) value from
% actual data Y and model data F. The code uses a general version of 
% R-square, based on comparing the variability of the estimation errors 
% with the variability of the original values. RSQUARE also outputs the
% root mean squared error (RMSE) for the user's convenience.
%
% Note: RSQUARE ignores comparisons involving NaN values.
% 
% INPUTS
%   Y       : Actual data
%   F       : Model fit
%
% OPTION
%   C       : Constant term in model
%             R-square may be a questionable measure of fit when no
%             constant term is included in the model.
%   [DEFAULT] TRUE : Use traditional R-square computation
%            FALSE : Uses alternate R-square computation for model
%                    without constant term [R2 = 1 - NORM(Y-F)/NORM(Y)]
%
% OUTPUT 
%   R2      : Coefficient of determination
%   RMSE    : Root mean squared error
%
% EXAMPLE
%   x = 0:0.1:10;
%   y = 2.*x + 1 + randn(size(x));
%   p = polyfit(x,y,1);
%   f = polyval(p,x);
%   [r2 rmse] = rsquare(y,f);
%   figure; plot(x,y,'b-');
%   hold on; plot(x,f,'r-');
%   title(strcat(['R2 = ' num2str(r2) '; RMSE = ' num2str(rmse)]))
%   
% Jered R Wells
% 11/17/11
% jered [dot] wells [at] duke [dot] edu
%
% v1.2 (02/14/2012)
%
% Thanks to John D'Errico for useful comments and insight which has helped
% to improve this code. His code POLYFITN was consulted in the inclusion of
% the C-option (REF. File ID: #34765).

if isempty(varargin); c = true; 
elseif length(varargin)>1; error 'Too many input arguments';
elseif ~islogical(varargin{1}); error 'C must be logical (TRUE||FALSE)'
else c = varargin{1}; 
end

% Compare inputs
if ~all(size(y)==size(f)); error 'Y and F must be the same size'; end

% Check for NaN
tmp = ~or(isnan(y),isnan(f));
y = y(tmp);
f = f(tmp);

if c; r2 = max(0,1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));
else r2 = 1 - sum((y(:)-f(:)).^2)/sum((y(:)).^2);
    if r2<0
    % http://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
        warning('Consider adding a constant term to your model') %#ok<WNTAG>
        r2 = 0;
    end
end

rmse = sqrt(mean((y(:) - f(:)).^2));
