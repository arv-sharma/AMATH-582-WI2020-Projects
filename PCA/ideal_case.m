%% LOADING DATA

clear; clc; close all
addpath('D:\course_work\amath582\hw3\data\')

C_1(3) = struct('data', [], 'results', []);
for jj = 1:1:3
    temp = load(strcat('cam', num2str(jj), '_1.mat'));
    temp1 = fieldnames(temp);
    C_1(jj).data = temp.(temp1{1});
end
clearvars temp

%% TRACKING THE PAINT CAN

for kk = 1:1:3
    
    sz = size(C_1(kk).data);
    C_1(kk).results.images = zeros(sz);
    C_1(kk).results.images = uint8(C_1(kk).results.images);
    x = 1:1:sz(2);
    y = 1:1:sz(1);
    num_of_frames = sz(end);
    [X, Y] = meshgrid(x, y);
    rec_mask = zeros(sz(1:2));
    if kk == 1
        rec_mask(X > 285 & X < 421) = 1;    % Restrict to interested area
    elseif kk == 2
        rec_mask(X > 220 & X < 383) = 1;
    else
        rec_mask(X > 221 & X < 469 & Y > 211 & Y < 335) = 1;
    end
    if kk ~= 3  
        circ_filt = fspecial('disk', 4);    % A circular disk filter
    else
        circ_filt = fspecial('disk', 2);
    end
    C_1(kk).results.x_loc = nan(sz(end), 1);
    C_1(kk).results.y_loc = nan(sz(end), 1);
    figure, imshow(C_1(kk).data(:, :, :, 1))
    % Get initial (x, y) pixel location from first frame
    prompt1 = 'Enter x-coordinate of torchlight from the frame: ';
    C_1(kk).results.x_loc(1) = input(prompt1);
    prompt2 = 'Enter y-coordinate of torchlight from the frame: ';
    C_1(kk).results.y_loc(1) = input(prompt2);
    close;

    for ii = 2: 1: num_of_frames

        I = C_1(kk).data(:, :, :, ii);
        I_gs = rgb2gray(I);
        I_gs_bl = imfilter(double(I_gs), circ_filt, 0);
        I_gs_bl_r = I_gs_bl .* rec_mask;    
        % dist_mask filters out pixels too far from previous known location
        dist_mask = sqrt((X - C_1(kk).results.x_loc(ii - 1)) .^ 2 +...
            (Y - C_1(kk).results.y_loc(ii - 1)) .^ 2);
        if kk ~= 3
            dist_mask(dist_mask < 19) = 1;
            dist_mask(dist_mask >= 19) = 0;
        else
            dist_mask(dist_mask < 30) = 1;
            dist_mask(dist_mask >= 30) = 0;
        end

        I_gs_bl_rd = I_gs_bl_r .* dist_mask;

        [~, maxpt_temp] = max(I_gs_bl_rd(:));

        I_gs_bl_rd_mod = I_gs_bl_rd;
        if kk ~= 3
            I_gs_bl_rd_mod(Y > (Y(maxpt_temp) - 8)) = 0;
        else
            I_gs_bl_rd_mod(X > (X(maxpt_temp) - 15)) = 0;
            I_gs_mod = double(I_gs) .* dist_mask .* rec_mask;
            I_gs_mod(X > (X(maxpt_temp) - 15)) = 0;
            [max_temp2, maxpt_temp2] = max(I_gs_mod(:));
        end

        [max_temp1, maxpt_temp1] = max(I_gs_bl_rd_mod(:));
        
        if max_temp1 > (.95 * 255)
            maxpt = maxpt_temp1;
        elseif kk == 3 && max_temp2 > (0.85 * 255)
            maxpt = maxpt_temp2;
        else
            maxpt = maxpt_temp;
        end    

        % Shading the detected point in red
        temp_img = I;
        temp_img((Y(maxpt) - 2):(Y(maxpt) + 2),...
            (X(maxpt) - 2):(X(maxpt) + 2), 1) = 255;
        temp_img((Y(maxpt) - 2):(Y(maxpt) + 2),...
            (X(maxpt) - 2):(X(maxpt) + 2), 2:3) = 0;
        C_1(kk).results.images(:, :, :, ii) = temp_img;

        [C_1(kk).results.x_loc(ii, 1), C_1(kk).results.y_loc(ii, 1)] =...
            deal(X(maxpt), Y(maxpt));

    end
    
end

%% PCA AND PLOTTING

[min_len, min_len_case] = min([length(C_1(1).results.x_loc),...
    length(C_1(2).results.x_loc), length(C_1(3).results.x_loc)]);
new_len = length(C_1(1).results.x_loc(10:end));

X_pca = [C_1(1).results.x_loc(10:end)';
         C_1(1).results.y_loc(10:end)';
         C_1(2).results.x_loc(19:(19 + new_len - 1))';
         C_1(2).results.y_loc(19:(19 + new_len - 1))';
         C_1(3).results.x_loc(9:(9 + new_len - 1))';
         C_1(3).results.y_loc(9:(9 + new_len - 1))'];

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
legend('Mode 1', 'Mode 2')
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
legend('Mode 1', 'Mode 2')
axis tight
xlabel({'Time, t', '(e)'})
ylabel({'Projection in', 'PCA mode'})
