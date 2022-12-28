function NoiseComposition

load('.\image\ValidationGtBlocksRaw.mat');
load('.\image\ValidationNoisyBlocksRaw.mat')
gt = reshape(ValidationGtBlocksRaw,[],256,256);
gt = shiftdim(gt,1);
noisy = reshape(ValidationNoisyBlocksRaw,[],256,256);
noisy = shiftdim(noisy,1);

for i = 1 : size(gt,3)
    rgb = demosaic(uint16(round(gt(:,:,i)*65535)),'grbg');
    rgb = single(rgb)/65535;
    luma = (rgb(:,:,1)+rgb(:,:,2)*2+rgb(:,:,3))/4;
    c = mean2((gt(:,:,i) - luma).^2);
    n = mean2((gt(:,:,i) - noisy(:,:,i)).^2);
    snr(i) = 10*log10(mean2(gt(:,:,i).^2)/n);
    r(i) = c/(c+n);
end
axes1 = axes('Parent',figure);
plot(snr,r,'.k');
ylabel('CNR');
xlabel('SNR(dB)');
set(axes1,'FontSize',16);
end