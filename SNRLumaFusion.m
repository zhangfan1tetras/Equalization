function SNRLumaFusion
h = [1 2 1;2 4 2;1 2 1]/16;
CCM = inv([ 1/4+1 1/2   1/4;...
            1/4   1/2+1 1/4;...
            1/4   1/2   1/4+1]/2);
load('.\image\ValidationGtBlocksRaw.mat');
load('.\image\ValidationNoisyBlocksRaw.mat');
gt = reshape(ValidationGtBlocksRaw,[],256,256);
gt = shiftdim(gt,1);
noisy = reshape(ValidationNoisyBlocksRaw,[],256,256);
noisy = shiftdim(noisy,1);

for  i = 1 : size(gt,3)
    i
    luman = imfilter(noisy(:,:,i),h);   % Estimate L
    rgb0 = single(demosaic(uint16(gt(:,:,i)*65535),'grbg'));
    rgb = single(demosaic(uint16(noisy(:,:,i)*65535),'grbg'));
    mse = mean((rgb(:) - rgb0(:)).^2);
    
    eql = (noisy(:,:,i) + luman)/2;     % Equalization at lamda = 2
    eql = double(demosaic(uint16(round(eql*65535)),'grbg'));
    eql = reshape(eql,[256*256 3]);
    eql = single(eql)*CCM';             % resaturation
    eql = reshape(eql,[256 256 3]);
    mse_eql = mean((eql(:) - rgb0(:)).^2);
    gain(i) = 10*log10(mse/mse_eql);
    snr(i) = 10*log10(mean2(gt(:,:,i).^2)/mean2((gt(:,:,i) - noisy(:,:,i)).^2));
end
axes1 = axes('Parent',figure);
plot(snr,gain,'.k');
ylabel('RGB-SNR gain(dB)');
xlabel('SNR(dB)');
set(axes1,'FontSize',16);
end
