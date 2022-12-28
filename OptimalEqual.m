function OptimalEqual
h = [1 2 1;2 4 2;1 2 1]/16;
CCM(:,:,1) = eye(3);
CCM(:,:,2) = inv([ 1/4+1 1/2   1/4;...
            1/4   1/2+1 1/4;...
            1/4   1/2   1/4+1]/2);
CCM(:,:,3) = inv([ 2/4+1 2/2   2/4;...
            2/4   2/2+1 2/4;...
            2/4   2/2   2/4+1]/3);
CCM(:,:,4) = inv([ 3/4+1 3/2   3/4;...
            3/4   3/2+1 3/4;...
            3/4   3/2   3/4+1]/4);
CCM(:,:,5) = inv([ 4/4+1 4/2   4/4;...
            4/4   4/2+1 4/4;...
            4/4   4/2   4/4+1]/5);
load('.\image\ValidationGtBlocksRaw.mat');
load('.\image\ValidationNoisyBlocksRaw.mat');
gt = reshape(ValidationGtBlocksRaw,[],256,256);
gt = shiftdim(gt,1);
noisy = reshape(ValidationNoisyBlocksRaw,[],256,256);
noisy = shiftdim(noisy,1);

for  i = 1 : size(gt,3)
    i
    luman = imfilter(noisy(:,:,i),h);
    rgb0 = single(demosaic(uint16(round(gt(:,:,i)*65535)),'grbg'));
    yuv0 = single(rgb2ycbcr(uint16(rgb0)));
    rgb = single(demosaic(uint16(round(noisy(:,:,i)*65535)),'grbg'));
    yuv = single(rgb2ycbcr(uint16(rgb)));
    mse = mean((rgb(:) - rgb0(:)).^2);
    mse_y = mean2((yuv(:,:,1) - yuv0(:,:,1)).^2);
    mse_uv = mean2((yuv(:,:,2)-yuv0(:,:,2)).^2)+mean2((yuv(:,:,3)-yuv0(:,:,3)).^2);
    snr(i) = 10*log10(mean2(gt(:,:,i).^2)/mean2((gt(:,:,i) - noisy(:,:,i)).^2));
    for j = 1 : 5
        eql = (noisy(:,:,i) + luman*(j-1))/j;
        eql = double(demosaic(uint16(eql*65535),'grbg'));
        eql = reshape(eql,[256*256 3]);
        eql = single(eql)*CCM(:,:,j)';
        eql = reshape(eql,[256 256 3]);
        yuv_eql = single(rgb2ycbcr(uint16(eql)));
        mse_eql = mean((eql(:) - rgb0(:)).^2);
        gain(i,j) = 10*log10(mse/mse_eql);
        gain_y(i,j) = 10*log10(mse_y/mean2((yuv_eql(:,:,1)-yuv0(:,:,1)).^2));
        gain_uv(i,j) = 10*log10(mse_uv/(mean2((yuv_eql(:,:,2)-yuv0(:,:,2)).^2)+mean2((yuv_eql(:,:,3)-yuv0(:,:,3)).^2)));
        gain_yuv(i,j) = 10*log10((mse_y+mse_uv)/(3*mean((yuv_eql(:)-yuv0(:)).^2)));
    end
    [max_gain(i),max_lamda(i)] = max(gain(i,:));
end
axes1 = axes('Parent',figure);
plot(snr,max_gain,'.k');
ylabel('RGB-SNR gain(dB)');
xlabel('SNR(dB)');
ylim(axes1,[-10 4]);
set(axes1,'FontSize',16);

axes2 = axes('Parent',figure);
plot(snr,max_lamda,'.k');
ylabel('Optimal \lambda');
xlabel('SNR(dB)');
ylim(axes2,[0 5]);
set(axes2,'FontSize',16);

axes3 = axes('Parent',figure);
plot3 = plot(1:5,mean(gain_yuv,1),'-hk',1:5,mean(gain_y,1),'-^k',1:5,mean(gain_uv,1),'-vk');
set(plot3(1),'DisplayName','YUV-SNR');
set(plot3(2),'DisplayName','     Y-SNR');
set(plot3(3),'DisplayName','  UV-SNR');
ylabel('SNR gain(dB)');
xlabel('\lambda');
set(axes3,'FontSize',16);
legend3 = legend(axes3,'show');
set(legend3,'Position',[0.15 0.20 0.24 0.18]);
end
