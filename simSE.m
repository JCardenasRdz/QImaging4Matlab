function Signal=simSE(pars,TR)
% Simulation of Spin Echo signal
%
% Author:
% Julio Cárdenas-Rodríguez
%  University of of Arizona
%
%  Sintax
%   Signal=simSE(pars,TR)
%   where pars=[Mz,T1];
%  
%   Mz          => Magnetization at thermal equilibrium
%   T1 `        => T1 time
%   TR          => Repetition time
%
R1=1/pars(2);
Signal=  pars(1) .* (1-exp(-R1.*TR));
end
