function [ sag_img, cor_img, hor_img ] = three_cut_maker(image_3d, num_cuts )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

sag_slices = (size(image_3d,1));
cor_slices = (size(image_3d,2));
hor_slices = size(image_3d,3);

sag_cuts = floor(sag_slices/num_cuts);
cor_cuts = floor(cor_slices/num_cuts);
hor_cuts = floor(hor_slices/num_cuts);

    sag_img = [];
    for i = 1:sag_cuts:sag_slices
    sag_img = horzcat(sag_img, rot90(squeeze(image_3d(i,:,:))));
    end
    
    cor_img = [];
    for j = 1:cor_cuts:cor_slices
    cor_img = horzcat(cor_img, rot90(squeeze(image_3d(:,j,:))));
    end
    
    hor_img = [];
    for k = 1:hor_cuts:hor_slices
    hor_img = horzcat(hor_img, rot90(squeeze(image_3d(:,:,k))));
    end

end

