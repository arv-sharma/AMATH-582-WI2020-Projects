clear
close all

%% Adding file path

addpath('D:\research_UW\PIT\AVI_Files\')
addpath(genpath('D:\research_UW\lib'))

%% Getting info about the file from log

prompt = '\n Enter filename: ';
filename = input(prompt, 's');
% filename = 'pit_0451';
temp = strsplit(filename, '_');
case_id = temp{2};
log_data = readtable('PIT_log.csv');
scopedata_filename = log_data.Shot_(strcmp(filename, log_data.KiranaFile));
scopedata_filename = scopedata_filename{:};
[~, case_idx] = whosi(log_data.KiranaFile, case_id);
Fs_img = log_data.KiranaFrameRate_fps_(case_idx);

%% Opening AVI and log files

[frames, num_of_frames, frame_size] = readavi(strcat(filename, '.avi'));
DPO_data = load(strcat('D:\research_UW\PIT\Scope_data\DPO_5054\mat\setup_test_',...
    case_id));
TDS_data = load(strcat('D:\research_UW\PIT\Scope_data\TDS_3034B\TDS_3034Bmat\setup_test_',...
    case_id));

%% Resampling data at each image time point

Fs_TDS = 1 / mode(diff(TDS_data.time));
% timeDelta_btw_img_in_TDS = Fs_TDS / Fs_img;
mod_TDS_timing_signal = TDS_data.CH1 - 1.5; % Because 1.5 V is where an image is triggered
% Next step: Use the diff function to find the points that are "rising",
% and use logical AND with a product of adjacent points. The product of
% adjacent points would be negative in this case.
img_pts_filt = find((diff(TDS_data.CH1) > 0) & ((mod_TDS_timing_signal(1:end - 1) .*...
    mod_TDS_timing_signal(2:end)) < 0));
if length(img_pts_filt) ~= 181
    error('Check the trigger signal.')
end
img_pts_TDS = img_pts_filt(1:180) + 1;

oldTDSfieldnames = fieldnames(TDS_data);
newTDSfieldnames = cellfun(@(x) strcat('img_', x), oldTDSfieldnames, 'UniformOutput', 0);
for ii = 1: 1: length(oldTDSfieldnames)
    TDS_data.(newTDSfieldnames{ii}) = TDS_data.(oldTDSfieldnames{ii})(img_pts_TDS);
end

oldDPOfieldnames = fieldnames(DPO_data);
newDPOfielnames = cellfun(@(x) strcat('img_', x), oldDPOfieldnames, 'UniformOutput', 0);

[~, ~, img_pts_DPO] = intersect(TDS_data.img_time, DPO_data.time);
for ii = 1: 1: length(oldDPOfieldnames)
    DPO_data.(newDPOfielnames{ii}) = DPO_data.(oldDPOfieldnames{ii})(img_pts_DPO);
end

Fs_DPO = 1 / mode(diff(DPO_data.time));

%% Convert to grayscale

for ii = 1:1:num_of_frames
    frames(ii).gray_pxdata = rgb2gray(frames(ii).pxdata);
    % Cropping
    frames(ii).cr_gray_pxdata = imcrop(frames(ii).gray_pxdata, [183, 116, 553, 553]);
end

%% Overall Intensity Analysis and Signal Analysis

intensity = nan(num_of_frames, 1);

for ii = 1:1:num_of_frames
    intensity(ii, 1) = sum(sum(frames(ii).gray_pxdata));
end

% Finding peaks

[pks1, locs1] = findpeaks(intensity, 'MinPeakHeight', 0.75 * max(intensity));
[pks2, locs2] = findpeaks(-intensity(1:locs1(1)));
start_of_int_phnmn = locs2(end);

% FFT

% Limages = num_of_frames - start_of_int_phnmn + 1;
Limages = num_of_frames;
% Fs_img = 5000000;
% Y = fft(intensity(start_of_int_phnmn:end));
Y = fft(intensity);
Timages = (0: 1 / Fs_img: (Limages - 1) / Fs_img)';
Freqimages = ((0: Limages-1) * (Fs_img / Limages))';
ImageFFT = abs(Y) .* (2 / Limages);
ImageFFT(1, 1) = ImageFFT(1, 1) / 2;

figure, plot(intensity)
hold on
plot(start_of_int_phnmn, -pks2(end), 'ro')

figure, plot(Freqimages, ImageFFT)

L_sig = length(DPO_data.time);
Y = fft(DPO_data.CH4);
freq_sig = ((0: L_sig - 1) * (Fs_DPO / L_sig));
SigFFT = abs(Y) .* (2 / L_sig);
SigFFT(1, 1) = SigFFT(1, 1) / 2;

figure
subplot(2, 1, 1)
plot(DPO_data.time, DPO_data.CH4)
subplot(2, 1, 2)
plot(freq_sig, SigFFT)
xlim([0 1e6])

%% SVD/PCA on images

Image_composite = nan(554 * 554, length([25:50]));
c1 = 0;
for frame_num = 25:75
    c1 = c1 + 1;
    frame = frames(frame_num).cr_gray_pxdata;
    Image_composite(:, c1) = frame(:);  
end

avg_img = mean(Image_composite, 2);
Image_composite_pca = Image_composite - avg_img * ones(1, length([25:50]));
[U, S, V] = svd(Image_composite_pca, 'econ');

figure, imagesc(reshape(avg_img, 554, 554))
figure;
for ii = 1: 1: 6
    subplot(2, 3, ii)
    imagesc(reshape(U(:, ii), 554, 554))
end
colormap(gray)

%% Getting a bright image for finding boundaries

cr_frame_size = size(frames(1).cr_gray_pxdata);
bright_sum_img = zeros(cr_frame_size);
for jj = locs1'
    bright_sum_img = bright_sum_img + double(frames(jj).cr_gray_pxdata);
end
bright_avg_img = uint8(bright_sum_img / length(locs1));    

%% Finding center of everything

% Assumption - The center is the same for all images in a video

fig = figure;
imshow(bright_avg_img)
I = bright_avg_img;
[centers_i, radii_i, metric_i] = imfindcircles(I, [40, 60], 'ObjectPolarity',...
    'dark', 'Sensitivity', 0.95);
if length(radii_i) == 1
    inner_radius = radii_i;
    image_center = centers_i;
    hold on
    plot(image_center(1), image_center(2), 'r*')
    viscircles(image_center, inner_radius, 'Color', 'm');
else
    error('Inner circle detection filed. Troubleshoot manually!')
end
[centers_o, radii_o, metric_o] = imfindcircles(I, [250, 280], 'ObjectPolarity',...
    'bright', 'Sensitivity', 0.98);
if length(radii_o) == 1
    outer_radius = radii_o;
    viscircles(image_center, outer_radius, 'Color', 'b');
else
    error('Outer circle detection filed. Troubleshoot manually!')
end

% Find if the edge detection is acceptable
figure(fig)
fprintf('\nDistance between centers of the two circles in image = %.1f pixels\n',...
    pdist([image_center; centers_o], 'euclidean'))
prompt1 = 'Enter 1 to continue, else abort: ';
if input(prompt1) ~= 1
    error('Program terminated.')
end

% T = adaptthresh(I, 'NeighborhoodSize', 4 * floor(size(I) / 20) + 1);
% % T = adaptthresh(I, 'NeighborhoodSize', 4 * floor(size(I) / 16) + 1);  %
% % for pit 0615
% I1 = imbinarize(I, T);
% figure, imshow(I1)
% [B, L] = bwboundaries(I1, 'noholes');
% boundL = cellfun(@length, B, 'UniformOutput', 0); % Lengths of extracted boundaries
% boundLvec = cell2mat(boundL);
% [~, boundLvecsort] = sort(boundLvec, 'descend');
% 
% % Outer boundary - For estimate of center
% for iii = 1: 1: length(boundLvecsort)
%     temp = B{boundLvecsort(iii)};
%     xmin = min(temp(:,2));
%     xmax = max(temp(:,2));
%     x0est = mean([xmin xmax]);
%     ymax = max(temp(:,1));
%     ymin = min(temp(:,1));
%     y0est = mean([ymin, ymax]);
%     if (abs(x0est - frame_size(2) / 2) < 0.05 * frame_size(2) / 2) && ...
%             (abs(y0est - frame_size(1) / 2) < 0.05 * frame_size(1) / 2) ...
%             && (xmin > 170 && xmin < 210)
%         ob_id = boundLvecsort(iii);
%         break
%     end
% end
%         
% outer_bound = B{ob_id};  % since the longest boundary is probably the outer boundary
% hold on
% plot(outer_bound(:, 2), outer_bound(:, 1), 'b-')
% 
% % ymax = max(outer_bound(:,1));
% % ymin = min(outer_bound(:,1));
% % y0est = mean([ymin, ymax]);
% % 
% % xmin = min(outer_bound(:,2));
% % xmax = max(outer_bound(:,2));
% % x0est = mean([xmin xmax]);
% 
% rsquare = sqrt((outer_bound(:,1) - y0est).^2 + (outer_bound(:,2) - x0est).^2);
% rest = max(rsquare);
% outer_pts_filt = rsquare > (0.98 * rest);
% 
% Outer_Bound = [outer_bound(outer_pts_filt, 1) outer_bound(outer_pts_filt,2)];
% hold on
% plot(Outer_Bound(:,2), Outer_Bound(:,1), 'r--')
% 
% % Fitting a non-linear least square fit to get center
% xdat = Outer_Bound(:,2);
% ydat = Outer_Bound(:,1);
% ytar = ydat.^2;
% lsfun = @(x) x(1).^2 - (xdat-x(2)).^2 + x(3).*(2*ydat - x(3)) - ytar;  
% From equation of a circle
% X0 = [rest, x0est, 0.7 * y0est];
% lb = [0.95*rest, .4*frame_size(2), .4*frame_size(1)];
% ub = [1.05*rest, .6*frame_size(2), .6*frame_size(1)];
% 
% X = lsqnonlin(lsfun, X0, lb, ub);
% outer_radius = X(1);
% % hold on
% % plot(X(2), X(3), 'r*')
% 
% % Inner boundary
% 
% for jjj = iii + 1: 1: length(boundLvecsort)
%     temp = B{boundLvecsort(jjj)};
%     ri = sqrt((temp(:,1) - y0est).^2 + (temp(:,2) - x0est).^2);
%     if all(ri < outer_radius)
%         ib_id = boundLvecsort(jjj);
%         break
%     end
% end
% % inner_bound = B{boundLvecsort(2)};  % since the 2nd longest boundary is
% probably the inner boundary
% inner_bound = B{ib_id};
% hold on
% plot(inner_bound(:, 2), inner_bound(:, 1), 'b-')
% % yimax = max(inner_bound(:,1));
% % yimin = min(inner_bound(:,1));
% % y0iest = mean([yimin, yimax]);
% % 
% % ximin = min(inner_bound(:,2));
% % ximax = max(inner_bound(:,2));
% % x0iest = mean([ximin ximax]);
% 
% % risquare = sqrt((inner_bound(:,1) - y0iest).^2 + (inner_bound(:,2) - x0iest).^2);
% riest = min(ri);
% inner_pts_filt = ri < (1.4 * riest);
% 
% Inner_Bound = [inner_bound(inner_pts_filt, 1) inner_bound(inner_pts_filt,2)];
% hold on
% plot(Inner_Bound(:,2), Inner_Bound(:,1), 'm-')
% 
% yimax = max(Inner_Bound(:,1));
% yimin = min(Inner_Bound(:,1));
% y0iest = mean([yimin, yimax]);
% 
% ximin = min(Inner_Bound(:,2));
% ximax = max(Inner_Bound(:,2));
% x0iest = mean([ximin ximax]);
% 
% % Fitting a non-linear least square fit to get center
% xidat = Inner_Bound(:,2);
% yidat = Inner_Bound(:,1);
% yitar = yidat.^2;
% lsfun = @(x) x(1).^2 - (xidat-x(2)).^2 + x(3).*(2*yidat - x(3)) - yitar;  
% From equation of a circle
% Xi0 = [riest, x0iest, 0.7 * y0iest];
% lb = [0.95*riest, .4*frame_size(2), .4*frame_size(1)];
% ub = [1.05*riest, .6*frame_size(2), .6*frame_size(1)];
% 
% Xi = lsqnonlin(lsfun, Xi0, lb, ub);
% hold on
% plot(Xi(2), Xi(3), 'mx')
% 
% inner_radius = Xi(1);
% image_center = Xi(2:3);

%% Sectorizing the images


I = frames(100).cr_gray_pxdata;
figure, imshow(I)
hold on
[X_px, Y_px] = meshgrid((1:1:cr_frame_size(2)), (1:1:cr_frame_size(1)));
X_diff_px = X_px - image_center(1);
Y_diff_px = -(Y_px - image_center(2));   % Flipping here because y = 1 is 
% actually the top row in image, not bottom
R_px = sqrt(X_diff_px .^ 2 + Y_diff_px .^ 2);
Theta_unadj_px = atan2(Y_diff_px, X_diff_px);
Theta_px = Theta_unadj_px;
Theta_px(Theta_unadj_px < 0) = Theta_px(Theta_unadj_px < 0) + (2 * pi);
Theta_deg_px = rad2deg(Theta_px);
contour(X_px', Y_px', Theta_deg_px', 24, 'ShowText', 'on')
figure, contourf(X_px', Y_px', R_px', 'ShowText', 'on')
figure, contourf(X_px', Y_px', Theta_deg_px', 12, 'ShowText', 'on')

% Sector half-width
sector_hw = 360 / 300;
theta_deg = 0 : 360 / 600 : 360 - 360 / 600;

% To see what's happening
% k2 = I;
figure;
for ii = 1:1:25
    k2 = bright_avg_img;
    k2(Theta_deg_px >= (theta_deg(ii) - sector_hw) &...
        Theta_deg_px < (theta_deg(ii) + sector_hw) & R_px > inner_radius &...
        R_px < 0.4 * outer_radius) = 255;
    imshow(k2)
    pause(0.1)
end 

%% Normalization of images by mean intensity in interested area

intensity_a = nan(num_of_frames, 1);

for ii = 1:1:num_of_frames
    intensity_a(ii, 1) = sum(sum(frames(ii).cr_gray_pxdata(R_px > inner_radius...
        & R_px < 0.4 * outer_radius)));
end

for ii = 1:1:num_of_frames
    frames(ii).norm_pxdata = frames(ii).cr_gray_pxdata * (mean(intensity_a) /...
        intensity_a(ii));
end

%% BinSums and Filtering

% Visualizing bins to rearrange BinSum vectors

Theta_deg_px_adj = Theta_deg_px - 90;
Theta_deg_px_adj(Theta_deg_px_adj < 0) = Theta_deg_px_adj(Theta_deg_px_adj < 0) + 360;

figure;
k2 = bright_avg_img;
k2(Theta_deg_px_adj >= (theta_deg(27) - sector_hw) & Theta_deg_px_adj <...
    (theta_deg(27) + sector_hw) & R_px > inner_radius & R_px < 0.4 * outer_radius) = 255;
k2(Theta_deg_px_adj >= (theta_deg(574) - sector_hw) & Theta_deg_px_adj <...
    (theta_deg(574) + sector_hw) & R_px > inner_radius & R_px < 0.4 * outer_radius) = 255;
imshow(k2)

BinSums = zeros(length(theta_deg), num_of_frames);
BinSums_full = zeros(length(theta_deg), num_of_frames); % These are bin sums of bigger sections

for mm = 1:1:num_of_frames
    temp_img = frames(mm).norm_pxdata;
    for nn = 1:1:length(theta_deg)
        BinSums(nn, mm) = sum(temp_img(Theta_deg_px_adj >= (theta_deg(nn) -...
            sector_hw) & Theta_deg_px_adj < (theta_deg(nn) +...
            sector_hw) & R_px > inner_radius & R_px < 0.4 * outer_radius));
        BinSums_full(nn, mm) = sum(temp_img(Theta_deg_px_adj >=...
            (theta_deg(nn) - sector_hw) & Theta_deg_px_adj < (theta_deg(nn)...
            + sector_hw) & R_px > inner_radius & R_px < outer_radius));
    end
end

% Noticing that bins 1 to 26 and 575 to 600 cover the light blocker
BinSums_tr = BinSums(27:574, :);
BinSums_full_tr = BinSums_full(27:574, :);

% Quick visualization of BinSums_r
figure
hold on
for kk = 30: 1: 35
    plot(BinSums_tr(:, kk))
end

BinSums_tr_t = fft(BinSums_tr);
BinSums_full_tr_t = fft(BinSums_full_tr);
N = size(BinSums_tr_t, 1);
k = (2 * pi / (N)) * [0:(N / 2 - 1) (-N / 2):-1];
ks = fftshift(k);

filt = (exp(- 9 * k .^2))';

BinSums_tr_tf = BinSums_tr_t .* repmat(filt, [1 180]);

BinSums_tr_f = ifft(BinSums_tr_tf);

% Quick visualization of BinSums_tr_f
colorvec = (linspace(0, 0.8, 6))' .* [1, 1, 1];
figure
hold on
for kk = 30: 1: 35
    plot(BinSums_tr_f(:, kk) / max(BinSums_tr_f(:, kk)), 'Color', colorvec(kk - 29, :))
end
ylabel('BinSums normalized by max')
xlabel('Azimuthal bin location')
title('BinSums')

% Quick visualization of BinSums_tr_f
colorvec = (linspace(0, 0.8, 6))' .* [1, 1, 1];
figure
hold on
for kk = 30: 1: 35
    plot(ks, fftshift(abs(BinSums_full_tr_t(:, kk))) /...
        max(abs(BinSums_full_tr_t(:, kk))), 'Color', colorvec(kk - 29, :))
end
xlim([0 3.14])

%% Speed

prompt2 = 'Enter frame range of interest: ';
frames_of_int = input(prompt2);  % Change this from case to case
top_pks = nan(10, length(frames_of_int));
top_locs = nan(10, length(frames_of_int));
spk_pks = nan(2, length(frames_of_int));
spk_locs = nan(2, length(frames_of_int));
for kk = frames_of_int
    [allpks, alllocs] = findpeaks(BinSums_tr_f(:, kk), 'MinPeakDistance', 30);
    [temppks, tempids] = sort(allpks, 'descend');
    top_pks(:, kk - (frames_of_int(1) - 1)) = temppks(1:10);
    top_locs(:, kk - (frames_of_int(1) - 1)) = alllocs(tempids(1:10));
end

% Locating sectors where max intensity occurs across frames

c1 = 0;
temp_top_pks = top_pks;
temp_top_locs = top_locs;
for nn = 1: 1: 10
    
    [~, top_ids_descend] = sort(temp_top_pks(:), 'descend');
    temp_sel_sector = temp_top_locs(top_ids_descend(1));
    temp_ids = find((temp_top_locs > temp_sel_sector - 9) &...
        (temp_top_locs < temp_sel_sector + 9));
    if length(temp_ids) == length(frames_of_int)
        c1 = c1 + 1;
        spk_locs(c1, :) = temp_top_locs(temp_ids);
        spk_pks(c1, :) = temp_top_pks(temp_ids);
        temp_top_pks(temp_ids) = [];
        temp_top_locs(temp_ids) = [];
    else
        temp_top_pks(temp_ids) = [];
        temp_top_locs(temp_ids) = [];
        continue
    end
end

% Calculating angular speed

ang_spd = deg2rad(diff(spk_locs, 1, 2) * 360 / 600) * Fs_img; % in rad/s
ang_freq = ang_spd / (2 * pi);

max_ang_freq = zeros(size(ang_spd, 1), 1);

for ii = 1: 1: (length(frames_of_int) - 1)
    for jj = ii + 1: 1: length(frames_of_int)
        temp = deg2rad((spk_locs(:, jj) - spk_locs(:, ii)) * 360 / 600) *...
            (Fs_img / (jj - ii)) / (2 * pi);
        max_ang_freq = max(max_ang_freq, abs(temp));
    end
end        

% Preparation for movie frames

spk_locs_theta_rd = deg2rad(theta_deg(spk_locs) + 90 + 26 * 0.6);
spk_rho = 0.4 * outer_radius;
% spk_rho = inner_radius + rescale(spk_pks, inner_radius, outer_radius);
spk_x = image_center(1) + spk_rho .* cos(spk_locs_theta_rd);
spk_y = image_center(2) - spk_rho .* sin(spk_locs_theta_rd);


%% Getting a video of the spoke movement
% colorvec = parula(4);
% c2 = 0;
% imagecount = length(frames_of_int);
% M(imagecount) = struct('cdata', [], 'colormap', []);
% writerObj = VideoWriter(strcat(filename, '.avi'));
% giffilename = strcat(filename, '.gif');
% writerObj.FrameRate = 2;
% open(writerObj);
% 
% for kk = frames_of_int
%     c2 = c2 + 1;
%     fig = figure;
%     fig.WindowState = 'maximized';
%     s1 = subplot(4, 1, 1:3);
%     imshow(frames(kk).cr_gray_pxdata)
%     hold on
%     plot(image_center(1), image_center(2), 'ro')
%     for mm = 1: 1: 3
%         plot([image_center(1) spk_x(mm, kk - (frames_of_int(1) - 1))],...
%             [image_center(2) spk_y(mm, kk - (frames_of_int(1) - 1))],...
%             'Color', colorvec(mm, :), 'LineWidth', 3)
%     end
%     
%     s2 = subplot(4, 1, 4);
%     s2.Box = 'on';
%     h1 = plot(DPO_data.img_time, DPO_data.img_CH4);
%     hold on
%     h2 = plot(DPO_data.img_time(kk), DPO_data.img_CH4(kk), 'r.');
%     h2.MarkerSize = 16;
%     
%     M(c2) = getframe(fig);
%     writeVideo(writerObj, M(c2));
%     
%     im = frame2im(M(c2));
%     [N1, map] = rgb2ind(im, 256);
%     if c2 == 1
%         imwrite(N1, map, giffilename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.5);
% 	else
% 		imwrite(N1, map, giffilename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
%     end
%     
%     close(fig)
% end
% close(writerObj)

%% Uniformity Evolution

% Binsums_tr_t are the FFT of Binsums_tr

samp_freq = 2 * pi / 600;
xdft = BinSums_full_tr_t(1:(N / 2 + 1), :);
psdx = (1 / (samp_freq * N)) * abs(xdft) .^ 2;
psdx(2:end - 1, :) = 2 * psdx(2:end - 1, :);
freq1 = 0: samp_freq / N: samp_freq - samp_freq / N;

band_p = bandpower(BinSums_full_tr);
DC_power_dens = psdx(1, :) ./ band_p;
norm_DC_power_dens = DC_power_dens / max(DC_power_dens);

figure;
subplot(2, 1, 1)
plot(DPO_data.img_time, norm_DC_power_dens)
xlim([0 DPO_data.img_time(end)])
subplot(2, 1, 2)
plot(DPO_data.img_time, DPO_data.img_CH4, 'LineWidth', 2)
hold on
plot(DPO_data.time, DPO_data.CH4, 'LineWidth', 0.5)
xlim([0 DPO_data.img_time(end)])

% Spectrogram

fig3 = figure(3);
pcolor(1:180, freq1, abs(BinSums_full_tr_t)), 
shading interp
ylim([0 samp_freq / 2])
xlabel('Frame')
ylabel('Frequency')
% set(gca, 'Ylim', [-50 50], 'Fontsize', [14]) 
colormap(pink)
fig3.Units = 'inches';
fig3.Position = [0 4 6.75 5.0625];
fig3.PaperUnits = 'inches';
fig3.PaperSize = [6.75 5.0625];
ax = gca;
ax.FontSize = 12;
