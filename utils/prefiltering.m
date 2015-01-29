
function [img_out] = prefiltering(img)

img_uint8(:,:,1) = uint8(img(:,:,1) * 255); 
img_uint8(:,:,2) = uint8(img(:,:,2) * 255);
img_uint8(:,:,3) = uint8(img(:,:,3) * 255);
%img_uint8(:,:,4) = uint8(img(:,:,4) * 255);

img_ao = imareaopen(img_uint8,100,8);
img_ao_ac = imareaclose(img_ao,100,8);

img_out(:,:,1)= medfilt2(mat2gray(img_ao_ac(:,:,1),[0 255]),[3 3]);
img_out(:,:,2)= medfilt2(mat2gray(img_ao_ac(:,:,2),[0 255]),[3 3]);
img_out(:,:,3)= medfilt2(mat2gray(img_ao_ac(:,:,3),[0 255]),[3 3]);
%img_out(:,:,4)= medfilt2(mat2gray(img_ao_ac(:,:,4),[0 255]),[3 3]);

end