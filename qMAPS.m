function [MAPS,MAPScfit]=qMAPS(Images_Matrix,xdata,mask)
% Process imaging data with an arbiritary function
%
% [MAPS,MAPScfit]=MAPSMap(Images_Matrix,xdata,mask)
%%  Inputs
%   Images_Matrix=    3D Matrix of MR Images of size (N,M,K)
%   xdata=            Vector of values for xdata in LSQcurvefit of size K
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

% ### Process Data ###
%  1) Reshape mask to vectorize and preallocate maps
[rows,cols,~]=size(Images_Matrix);
mask=reshape(mask,rows*cols,[]) ;
MAPS.map=nan(size(mask));
MAPS.rsq=MAPS.map;

%  2)  Find indexes for pixeles to be analyzed
indices= find(mask);

%  3)  preallocate MAPScfit cell
MAPScfit=cell(size(indices));

%  5)  Reshape T1 vxdata 3D matrix into a 2D matrix
Images_Matrix=reshape(Images_Matrix,rows*cols,[]) ;

%  6) CurveFitting for each voxel
x0=[1.2,3.0];
lb=[1,1];
ub=[2,5];
N=length(xdata);

for q=1:length(indices);
    
    Signal= Images_Matrix (indices(q),:)';   Signal=Signal./max(Signal);

[beta,~,resid,~,~,~,J]=lsqcurvefit(MAPSfunc,...
                                              x0,xdata,Signal,lb,ub,options);
    
    T1    = beta(2);
    T1ci  = nlparci(beta,resid,'jacobian',J);
    rmse  =rsquare(Signal,  MAPSfunc(beta,xdata)   );
    
    MAPScfit{q}.cfit.T1=     T1;
    MAPScfit{q}.cfit.T1ci=   T1ci;
    MAPScfit{q}.cfit.MAPS.map=   MAPS.map;
    
    MAPS.rsq(indices(q),1)=rmse;
    MAPS.map(indices(q),1)=T1;
end

%  Reshapes MAPS back to square matrices

MAPS.map=reshape(MAPS.map,rows,cols); 
MAPS.rsq=reshape(MAPS.rsq,rows,cols);
