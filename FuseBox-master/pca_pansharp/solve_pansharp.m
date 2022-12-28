function [ sharpened ] = solve_pansharp( image, pan )

if (size(image, 1) ~= size(pan, 1))
    image = upsample_ms(image);
else
    image = demosaic(uint16(image*65535),'gbrg');
    image = single(image)/65535;
end 

[m, n, d] = size(image);

M = reshape(image, m*n, d);

%pca transform on ms bands
[PCAData, PCAMap] = multi_pca(M, d);
PCAData = reshape(PCAData, [m, n, d]);
F = PCAData;

P = (pan - mean(pan(:))) * std(F(:))/std(pan(:)) + mean(F(:));

%replace 1st band with pan
F(:,:,1) = P;

%inverse PCA
F = inv_pca(PCAMap.M, reshape(F, m*n, d), PCAMap.mean);

sharpened = reshape(F, [m, n, d]);
