function [ sharpened ] = solve_pansharp( image, pan )

if (size(image, 1) ~= size(pan, 1))
    image = upsample_ms(image);
else
    image = demosaic(uint16(image*65535),'gbrg');
    image = single(image)/65535;
end 

I = image(:,:,1)*1/3 + image(:,:,2)*1/3 + image(:,:,3)*1/3;

P = (pan - mean(pan(:))) * std(I(:))/std(pan(:)) + mean(I(:));

sharpened = zeros(size(image));

for n=1:3
    sharpened(:,:,n)= image(:,:,n) + P - I;
end

end

