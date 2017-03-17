function [tsnr_image] = tsnr_creator(image_name)
%tsnr_creator takes an input image and returns tsnr 
%   This function will take an input image string ending with nii/nii.gz
%   return an image of the tsnr, that is, mean / standard deviation

img_data = load_nii(image_name);
img_data = img_data.img;
mean_img = squeeze(mean(img_data,4));
std_img = squeeze(std(img_data,0,4));
tsnr_image = mean_img./std_img;

end

