% 
clear all; clc
load T1_vTR_example.mat
Images_Matrix=VTR(:,:,2:3:end);
%%
[ ROI] = noise_mask( Images_Matrix(:,:,end) );
%% Allocate MEthods Structure
Methods.function=@(xdata,Signal) fit_T1vtr(xdata,Signal);
Methods.mask=ROI;
Methods.NumMaps=5;
xdata=Repetition_Time';

tic
[MAPS]=sqMAPS(Images_Matrix,xdata,Methods);
toc
%%

 
