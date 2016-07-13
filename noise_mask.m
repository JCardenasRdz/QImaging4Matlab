function [ mask,probabilities] = noise_mask( Images )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[rows,cols,Zdim]=size(Images);

Number_of_clusters=2;
Images_Matrix=reshape(Images,[],Zdim);

    
    [~,U,~] = fcm(Images_Matrix,Number_of_clusters,[nan nan 200 0]);
    index_tissue=find(  U(:,1) ~= max(U(:,1)));
    
probabilities=reshape(U(index_tissue,:),rows,cols);
mask=probabilities>=0.5;
end

