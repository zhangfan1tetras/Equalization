function BayerwhiteFusion
m = 3472;
n = 4640;
bl = 64;
folder = '.\Samsung\20220625231441_iso_12231_1600_shutter_160000000_lux_0.125000\';
%folder = '.\Samsung\20220625232441_iso_12231_1600_shutter_160000000_lux_0.125000\';
fid = fopen([folder 'noisy_2.raw']);
I  = fread(fid,'uint16');
fclose(fid);
I = reshape(I,[2*n m]);
bayer = double(I(1:2:end,:)');
white = double(I(2:2:end,:)');
% imwrite(uint8((bayer-bl)*2),'bayerchannel.jpg','quality',100);
% imwrite(uint8((white-bl)*2),'whitechannel.jpg','quality',100);

fid = fopen([folder 'gt.raw']);
I  = fread(fid,'uint16');
fclose(fid);
I = reshape(I,[2*n m]);
bayer0 = single(I(1:2:end,:)');
white0 = single(I(2:2:end,:)');
clear I;

bayer1 = bayer0 - bl;
gmean =  mean2(bayer1(1:2:end,1:2:end)+bayer1(2:2:end,2:2:end))/2;
bayer1(1:2:end,2:2:end) = gmean/mean2(bayer1(1:2:end,2:2:end))*bayer1(1:2:end,2:2:end);
bayer1(2:2:end,1:2:end) = gmean/mean2(bayer1(2:2:end,1:2:end))*bayer1(2:2:end,1:2:end);
rgb0 = single(demosaic(uint16(bayer1),'gbrg'));

bayer1 = bayer - bl;
gmean =  mean2(bayer1(1:2:end,1:2:end)+bayer1(2:2:end,2:2:end))/2;
bayer1(1:2:end,2:2:end) = gmean/mean2(bayer1(1:2:end,2:2:end))*bayer1(1:2:end,2:2:end);
bayer1(2:2:end,1:2:end) = gmean/mean2(bayer1(2:2:end,1:2:end))*bayer1(2:2:end,1:2:end);
rgb = single(demosaic(uint16(bayer1),'gbrg'));
clear bayer1;
mse = mean((rgb(:) - rgb0(:)).^2);
imwrite(uint8(rgb0*2),'gt.jpg','quality',100);
imwrite(uint8(rgb*2),'noisy.jpg','quality',100);
%%
calib = cat(3,bayer0(2:2:end,1:2:end),(bayer0(1:2:end,1:2:end)+bayer0(2:2:end,2:2:end))/2,bayer0(1:2:end,2:2:end));
calib = reshape(calib-bl,[],3);
w0 = imresize(white0-bl,0.5);
M = inv(calib'*calib)*calib'*w0(:);
clear calib w0;
%%
yuv0 = single(rgb2ycbcr(uint16((rgb0/(1023-bl)).^1*(1023-bl))));
yuv = single(rgb2ycbcr(uint16((rgb/(1023-bl)).^1*(1023-bl))));
mse_y = mean2((yuv(:,:,1) - yuv0(:,:,1)).^2);
mse_uv = mean2((yuv(:,:,2)-yuv0(:,:,2)).^2)+mean2((yuv(:,:,3)-yuv0(:,:,3)).^2);
snr = 10*log10(mean((rgb0(:)-bl).^2)/mse);

%%
eql = (bayer + white)/2; % lamda = 2
eql = single(demosaic(uint16(round(eql)),'gbrg')-bl);
%imwrite(uint8(eql*2),'desaturation.jpg','quality',100);
lamda = 1.6;%2;%
a2 = lamda - 1;
CCM = inv([ a2*M(1)+1 a2*M(2)   a2*M(3);...
            a2*M(1)   a2*M(2)+1 a2*M(3);...
            a2*M(1)   a2*M(2)   a2*M(3)+1]/(1+a2));

eql = reshape(eql,[m*n 3]);
eql = single(eql)*CCM';
eql = reshape(eql,[m n 3]);
eql(:,:,1) = mean2(eql(:,:,2))/mean2(eql(:,:,1))*eql(:,:,1);
eql(:,:,3) = mean2(eql(:,:,2))/mean2(eql(:,:,3))*eql(:,:,3);
eql = max(eql,0);
eql = mean(rgb0(:)) - mean(eql(:)) + eql;
imwrite(uint8(eql*2),['equalization' num2str(a2) '.jpg'],'quality',100);

yuv_eql = single(rgb2ycbcr(uint16(round((max(eql,0)/(1023-bl)).^1*(1023-bl)))));
mse_eql = mean((eql(:) - rgb0(:)).^2);
gain = 10*log10(mse/mse_eql);
gain_y = 10*log10(mse_y/mean2((yuv_eql(:,:,1)-yuv0(:,:,1)).^2));
gain_uv = 10*log10(mse_uv/(mean2((yuv_eql(:,:,2)-yuv0(:,:,2)).^2)+mean2((yuv_eql(:,:,3)-yuv0(:,:,3)).^2)));
gain_yuv = 10*log10((mse_y+mse_uv)/(3*mean((yuv_eql(:)-yuv0(:)).^2)));


%% Pansharpen
%bayer0 = bayer0((m/4+1):3*m/4,(n/4+1):3*n/4);
%white0 = white0((m/4+1):3*m/4,(n/4+1):3*n/4);
%bayer = bayer((m/4+1):3*m/4,(n/4+1):3*n/4);
%white = white((m/4+1):3*m/4,(n/4+1):3*n/4);
method = {'ground truth','Bayer only','equalization',...
    'guided','ihs','pca','wavelet','brovey','pxs'};
    %'brovey','wavelet','pca','guided','ihs','pxs',...
    %'closed','mmp'};% -7.7107 -4.4553 -2.0923 0.1018 -0.3681 -9.0399 0.1288 -1.0523...
% for i = 4 : length(method)
%     method{i}
%     sharpened = pansharpen(bayer,white,[method{i} '_panshar']);
%     sharpened = sharpened*1023 - bl; %???
%     sharpened(:,:,1) = mean2(sharpened(:,:,2))/mean2(sharpened(:,:,1))*sharpened(:,:,1);
%     sharpened(:,:,3) = mean2(sharpened(:,:,2))/mean2(sharpened(:,:,3))*sharpened(:,:,3);
%     sharpened = mean2(rgb0) - mean2(sharpened) + sharpened;
%     mse_shp = mean((sharpened(:) - rgb0(:)).^2);
%     gain(i) = 10*log10(mse/mse_shp)
%     sharpened = uint8((sharpened*2));
%     imwrite(sharpened,[method{i} '.jpg'],'quality',100);
% end
    
loc = [870 1150;1290 1530;1150 2210;1750 2210;2200 2050;500 1580]';    
% loc = [1150 2200  870 1750 1290  500;...
%        2210 2050 1150 2210 1530 1580];
siz = 256;
bdv = ones(siz+1,10,3)*255;
bdh = ones(10,(siz+1+10)*size(loc,2)+siz,3)*255;
I = [];
for i = 1 : length(method)
    ref = imread('gt.jpg');
    sharpened = imread([method{i} '.jpg']);
    row = ones(siz+1,siz,3)*255;
    for j = 1 : size(loc,2)
        patch = sharpened(loc(1,j):loc(1,j)+siz,loc(2,j):loc(2,j)+siz,:);
        row = [row bdv patch];
        refpt = single(ref(loc(1,j):loc(1,j)+siz,loc(2,j):loc(2,j)+siz,:));
        snr(i,j) = 10*log10(mean(refpt(:).^2)/mean((refpt(:)-single(patch(:))).^2));
    end
    I = [I;row;bdh];
end

figure1 = figure;
imshow(I,'Border','tight');
for i = 1 : length(method)
    annotation(figure1,'textbox',[0 1.03-0.111*i 0.16 0.04],...
        'String',{method{i}},'FontSize',8,'EdgeColor','none');
end
for i = 2 : length(method)
    for j = 1 : size(loc,2)
        annotation(figure1,'textbox',[0.05+0.14*j 0.999-0.111*i 0.0775 0.0341], ...
            'String',{[sprintf('%.1f',snr(i,j)) 'dB']},...
            'FontSize',6,'EdgeColor','none','Color',[1 1 1]);
    end
end


% I1 = imread('desaturation0.6.jpg');
% I2 = imread('desaturation1.jpg');
% I = [repmat([uint8((bayer-bl)*2) uint8((white-bl)*2)],[1 1 3]);...
%     I1 I2];
%imwrite(I,'result4.jpg','quality',100);
I = imread('result4.jpg');
figure2 = figure;
imshow(I,'Border','tight');
annotation(figure2,'textbox',[0.42 0.95 0.16 0.04],'Color',[1 1 1],...
        'String',{'Bayer'},'FontSize',10,'EdgeColor','none');
annotation(figure2,'textbox',[0.92 0.95 0.16 0.04],'Color',[1 1 1],...
        'String',{'white'},'FontSize',10,'EdgeColor','none');
annotation(figure2,'textbox',[0.33 0.45 0.18 0.04],'Color',[1 1 1],...
        'String',{'partial desaturation'},'FontSize',10,'EdgeColor','none');
annotation(figure2,'textbox',[0.85 0.45 0.16 0.04],'Color',[1 1 1],...
        'String',{'full desaturation'},'FontSize',10,'EdgeColor','none');

% pnshp = hcm(bayer,white)-bl;
% pnshp(:,:,1) = mean2(pnshp(:,:,2))/mean2(pnshp(:,:,1))*pnshp(:,:,1);
% pnshp(:,:,3) = mean2(pnshp(:,:,2))/mean2(pnshp(:,:,3))*pnshp(:,:,3);
% imwrite(uint8((pnshp)*2),['RGBW_pansharp.01.jpg'],'quality',100);
end

function sharpened = pansharpen(bayer,white,method)
addpath(['.\FuseBox-master\common\']);
addpath(['.\FuseBox-master\' method]);
%rgb = cat(3,bayer(2:2:end,1:2:end),(bayer(1:2:end,1:2:end)+bayer(2:2:end,2:2:end))/2,bayer(1:2:end,2:2:end));
white = white/1023;
bayer = bayer/1023;
sharpened = solve_pansharp(bayer, white);
end

function pnshp = hcm(bayer,white)
addpath('.\HCM package v2\');
HR_image = white;
LR_image = single(demosaic(uint16(bayer),'gbrg'));
patch_size = 2;
reg_param = 0.1;
pnshp = HCM(HR_image,LR_image,[],patch_size,reg_param);
end