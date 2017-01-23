function [MAPS]=sqMAPS(Images_Matrix,xdata,Methods)
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


% Extract information from Methods Structure
mask            =   Methods.mask;
Number_of_Maps  =   Methods.NumMaps;
Processing_Function=Methods.function;

% ### Process Data ###
%  1) Reshape mask to vectorize and preallocate maps
[rows,cols,~]=size(Images_Matrix);
mask=reshape(mask,rows*cols,[]) ;
MAPS=cell(Number_of_Maps);
for r=1:Number_of_Maps
MAPS{r}=nan(rows,cols);
end

%  2)  Find indexes for pixeles to be analyzed
indices= find(mask);

%  3)  Reshape T1 vxdata 3D matrix into a 2D matrix
Images_Matrix=reshape(Images_Matrix,rows*cols,[]) ;

%  4) Process each voxel
for q=1:length(indices);
    Signal= Images_Matrix (indices(q),:)';
    Pars=Processing_Function(xdata,Signal);
    disp(['Current voxel = ',num2str(q)])
    for r=1:Number_of_Maps
       MAPS{r}(indices(q)) = Pars(r);
    end    
end
