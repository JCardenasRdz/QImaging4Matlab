function [FittingPars] =  fit_T1vtr(RepetionTime,Signal)
% Estimated T1 relaxation time using variable TR data
%
%% Syntaxis
%
%   [FittingPars,Ypred,Yci,vtrModel] =  fit_T1_vtr(RepetionTime,Signal)
%
%% Inputs
%   RepetionTime= N x 1 vector of repetition time for VTR experiment. Units= seconds
%   Signal= N x 1 Matrix of signal for the VTR experiment. Units= AU
%
%% Ouputs
%   FittingPars= 5 X 1 numeric array with the following elements
%       FittingPars(1)= T1 time (estimated)
%       FittingPars(2)= standard eror of T1 time
%       FittingPars(3)= Lower 95% confience interval for the estimated T1
%       FittingPars(4)= Upper 95% confience interval for the estimated T1
%       FittingPars(5)= Adjusted Rsquared for the model 
%
%   Ypred= N X 1 numeric array of the predicted signal
%
%   Yci= N X 2 numeric array with the following elements
%       Yci(:,1)= Lower 95% confience interval for the estimated signal
%       Yci(:,2)= Upper 95% confience interval for the estimated signal
%
%% Example
%
%   TR=0:.1:6;      
%   RepetionTime=TR';         
%   T1=2.7;
%   S=100 * ( 1 -exp (-TR./T1))';   Signal=awgn(S,25,'measured');
%   [FittingPars,Ypred,Yci,vtrModel] =  fitT1vtr(RepetionTime,Signal);
%  
%     plot(RepetionTime,Signal,'or',RepetionTime,Ypred,'-',RepetionTime,Yci,'--k');
%     xlabel('Repetition Time (sec)');    ylabel('Signal (au)');
%     title(['Estimated T1= ',num2str(FittingPars(1)), '   R-squared= ', num2str(FittingPars(end))])
%     legend({'Data','Predicted','Upper Bound','Lower Bound'},'Location','Best')
%
%%  Author
% Julio Cárdenas-Rodríguez
% University of Arizona
% Tucson, AZ
% cardenaslab.org

% Normalize
Signal=Signal./max(Signal);

% Initial Guess
x0=[max(Signal)*1.2,2.0];

% function
SE=@(p,xdata)  p(1) .* (1-exp(-xdata/p(2)));


%
options = optimoptions('lsqcurvefit');
options.Display='off';
 FittingPars=lsqcurvefit(SE,x0,RepetionTime, Signal,[0,0],[],options);

% Fit variable TR model
vtrModel=fitnlm(RepetionTime,Signal,SE,FittingPars(1:2));

% Allocate Variables
T1pred=vtrModel.Coefficients.Estimate(2);
T1stde=vtrModel.Coefficients.SE(2);
ci = coefCI(vtrModel);
T1ci=ci(2,:);
Rsquared=vtrModel.Rsquared.Adjusted;
FittingPars=[T1pred,T1stde,T1ci,Rsquared]'; 
end

