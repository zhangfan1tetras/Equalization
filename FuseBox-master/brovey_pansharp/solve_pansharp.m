function [ sharpened ] = solve_pansharp( image, pan )
if size(image,1) ~= size(pan,1)
    image = upsample_ms(image);
else
    image = demosaic(uint16(image*65535),'gbrg');
    image = single(image)/65535;
end    
    

[~, ~, d_im] = size(image);

sharpened = zeros(size(image));

denom = sum(image, 3);

for k = 1 : d_im
   sharpened(:,:,k) = (d_im * pan .* image(:,:,k) ) ./ denom;
end

end