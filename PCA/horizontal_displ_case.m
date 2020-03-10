%% LOADING DATA

clear; clc; close all
addpath('D:\course_work\amath582\hw3\data\')

C_3(3) = struct('data', [], 'results', []);
for jj = 1:1:3
    temp = load(strcat('cam', num2str(jj), '_3.mat'));
    temp1 = fieldnames(temp);
    C_3(jj).data = temp.(temp1{1});
end
clearvars temp

%% TRACKING THE PAINT CAN

for kk = 1:1:3
    
    sz = size(C_3(kk).data);
    C_3(kk).results.images = zeros(sz);
    C_3(kk).results.images = uint8(C_3(kk).results.images);
    x = 1:1:sz(2);
    y = 1:1:sz(1);
    num_of_frames = sz(end);
    [X, Y] = meshgrid(x, y);
    rec_mask = zeros(sz(1:2));
    if kk == 1
        rec_mask(X > 270 & X < 380) = 1;    % Restrict to interested area
    elseif kk == 2
        rec_mask(X > 220 & X < 390) = 1;
    else
        rec_mask(X > 260 & X < 425 & Y > 185 & Y < 320) = 1;
    end
%     if kk ~= 3
%         circ_filt = fspecial('disk', 3);
    if kk == 3
        circ_filt = fspecial('disk', 2);    % A circular disk filter
    end
    C_3(kk).results.x_loc = nan(sz(end), 1);
    C_3(kk).results.y_loc = nan(sz(end), 1);
    figure, imshow(C_3(kk).data(:, :, :, 1))
    % Get initial (x, y) pixel location from first frame
    prompt1 = 'Enter x-coordinate of torchlight from the frame: ';
    C_3(kk).results.x_loc(1) = input(prompt1);
    prompt2 = 'Enter y-coordinate of torchlight from the frame: ';
    C_3(kk).results.y_loc(1) = input(prompt2);
    close;

    for ii = 2: 1: num_of_frames

        I = C_3(kk).data(:, :, :, ii);
        if kk ~= 3
            I_r = double(I) .* repmat(rec_mask, [1, 1, 3]);
            r2g = I_r(:, :, 1) ./ I_r(:, :, 2);
            r2b = I_r(:, :, 1) ./ I_r(:, :, 3);
            % pink_filt looks for pink colored pixels
            pink_filt = (r2g > 1.2 & r2g < 1.55 & r2b > 1.1 & r2b < 1.5);
        else
            I_gs = rgb2gray(I);
            I_gs_bl = imfilter(double(I_gs), circ_filt, 0);
            I_gs_bl_r = I_gs_bl .* rec_mask;
        end
        % dist_mask filters out pixels too far from previous known location
        dist_mask = sqrt((X - C_3(kk).results.x_loc(ii - 1)) .^ 2 +...
            (Y - C_3(kk).results.y_loc(ii - 1)) .^ 2);
        if kk == 1
            dist_mask(dist_mask < 20) = 1;
            dist_mask(dist_mask >= 20) = 0;
        elseif kk == 2
            dist_mask(dist_mask < 30) = 1;
            dist_mask(dist_mask >= 30) = 0;
        else
            dist_mask(dist_mask < 20) = 1;
            dist_mask(dist_mask >= 20) = 0;
        end
        
        if kk ~= 3
            pink_filt_d = pink_filt .* dist_mask;
            x_pt = round(sum(sum(pink_filt_d .* X)) / sum(pink_filt_d(:)));
            y_pt = round(sum(sum(pink_filt_d .* Y)) / sum(pink_filt_d(:)));
        else
            I_gs_bl_rd = I_gs_bl_r .* dist_mask;
            [~, maxpt_temp] = max(I_gs_bl_rd(:));

            I_gs_bl_rd_mod = I_gs_bl_rd;
            I_gs_bl_rd_mod(X > (X(maxpt_temp) - 15)) = 0;
            I_gs_mod = double(I_gs) .* dist_mask .* rec_mask;
            I_gs_mod(X > (X(maxpt_temp) - 15)) = 0;
            
            [max_temp1, maxpt_temp1] = max(I_gs_bl_rd_mod(:));
            [max_temp2, maxpt_temp2] = max(I_gs_mod(:));


            if max_temp2 > (0.97 * 255)
                maxpt = maxpt_temp2;
            elseif max_temp1 > (.95 * 255)
                maxpt = maxpt_temp1;
            else
                maxpt = maxpt_temp;
            end
            [x_pt, y_pt] = deal(X(maxpt), Y(maxpt));
        end

        % Shading the detected point in red
        temp_img = I;
        temp_img((y_pt - 2):(y_pt + 2),...
            (x_pt - 2):(x_pt + 2), 1) = 255;
        temp_img((y_pt - 2):(y_pt + 2),...
            (x_pt - 2):(x_pt + 2), 2:3) = 0;
        C_3(kk).results.images(:, :, :, ii) = temp_img;
        
        [C_3(kk).results.x_loc(ii, 1), C_3(kk).results.y_loc(ii, 1)] =...
            deal(x_pt, y_pt);
        
    end
    
end

%% MONTAGE

multi = cat(4, C_3(2).results.images(:, :, :, 2:18:56));
fig1 = figure;
montage(multi)
fig1.Units = 'inches';
fig1.Position = [-.1 1.8 6.75 5.0625];
fig1.PaperUnits = 'inches';
fig1.PaperSize = [6.75 5.0625];

%% PCA AND PLOTTING

[min_len, min_len_case] = min([length(C_3(1).results.x_loc),...
    length(C_3(2).results.x_loc), length(C_3(3).results.x_loc)]);
% new_len = length(C_3(3).results.x_loc(13:end));
new_len = 223;

X_pca = [C_3(1).results.x_loc(17:(17 + new_len - 1))';
         C_3(1).results.y_loc(17:(17 + new_len - 1))';
         C_3(2).results.x_loc(4:(4 + new_len - 1))';
         C_3(2).results.y_loc(4:(4 + new_len - 1))';
         C_3(3).results.x_loc(13:(13 + new_len - 1))';
         C_3(3).results.y_loc(13:(13 + new_len - 1))'];

X_pca_mod = X_pca - mean(X_pca, 2); % Standardizing
[U, S, V] = svd(X_pca_mod' / sqrt(new_len - 1));
sig = diag(S);
Y_pca = V' * X_pca_mod;

fig = figure;
fig.Units = 'inches';
fig.Position = [-.1 1.8 6.75 6];
fig.PaperUnits = 'inches';
fig.PaperSize = [6.75 6];

s1 = subplot(3, 2, 1);
s1.Box = 'on';
h1 = plot(sig, 'ko-', 'LineWidth', 1.1);
xlabel({'Mode, k'; '(a)'})
ylabel('Singular value, \sigma_k')
axis tight

s2 = subplot(3, 2, 2);
s2.Box = 'on';
h2 = plot(cumsum(sig) / sum(sig), 'k^-', 'LineWidth', 1.1);
xlabel({'Mode, k'; '(b)'})
ylabel({'Cumulative energy in', 'first k modes, \sigma_k'})
axis tight

s3 = subplot(3, 2, 3:4);
s3.Box = 'on';
hold on
h3 = plot(U(:, 1), 'k', 'LineWidth', 1.5);
h4 = plot(U(:, 2), 'k-.', 'LineWidth', 1.5);
h4a = plot(U(:, 3), 'k--', 'LineWidth', 1.5);
h4b = plot(U(:, 4), 'k:', 'LineWidth', 1.5);
legend('Mode 1', 'Mode 2', 'Mode 3', 'Mode4')
axis tight
xlabel({'Time, t', '(c)'})
ylabel('PCA mode variation')

s4 = subplot(3, 2, 5);
s4.Box = 'on';
hold on
h5 = plot(X_pca_mod(1, :), 'k', 'LineWidth', 1.5);
h6 = plot(X_pca_mod(2, :), 'k-.', 'LineWidth', 1.5);
h7 = plot(X_pca_mod(3, :), 'k--', 'LineWidth', 1.5);
h8 = plot(X_pca_mod(4, :), 'k', 'LineWidth', 2);
h9 = plot(X_pca_mod(5, :), 'k-.', 'LineWidth', 2);
h10 = plot(X_pca_mod(6, :), 'k--', 'LineWidth', 2);
legend('Cam1x', 'Cam1y', 'Cam2x', 'Cam2y', 'Cam3x', 'Cam3y')
axis tight
xlabel({'Time, t', '(d)'})
ylabel('Position in pixel co-ordinates')

s5 = subplot(3, 2, 6);
s5.Box = 'on';
hold on
h11 = plot(Y_pca(1, :), 'k', 'LineWidth', 1.5);
h12 = plot(Y_pca(2, :), 'k-.', 'LineWidth', 1.5);
h12a = plot(Y_pca(3, :), 'k--', 'LineWidth', 1.5);
h12b = plot(Y_pca(4, :), 'k:', 'LineWidth', 1.5);
legend('Mode 1', 'Mode 2', 'Mode3', 'Mode 4')
axis tight
xlabel({'Time, t', '(e)'})
ylabel({'Projection in', 'PCA mode'})
