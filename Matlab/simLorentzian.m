function Lorentzian=simLorentzian(pars,ppm)
% Simulation of Lorentzian function for fitting of CEST data
% Reference
% Saizz et a; Journal of Magnetic Resonance 211 (2011) 149?155
% Author:
% Julio Cárdenas-Rodríguez
%  University of of Arizona
%  Sintax
%   [L]=simLorentzian(pars,ppm)
%  where pars=[Amplitude,Width,Center];
%  Input
%  Properties of the Lorentzian function
%   Amplitude   => amplitude
%   Width `     => Width at half max
%   Center      => Offset at which the function is centered
%
%
Amplitude=pars(1);
        Gamma= pars(2)^2 /4;
            OffsetSqr= (ppm-pars(3)).^2;


Lorentzian= (Amplitude.*Gamma) ./ (Gamma+OffsetSqr);


end
