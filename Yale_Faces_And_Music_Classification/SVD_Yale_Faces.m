%% YALE FACES

%% LOAD CROPPED DATA
clear; clc; close all

addpath('D:\course_work\amath582\hw4\yalefaces_cropped\CroppedYale')
yale_cr_folders = dir(['D:\course_work\amath582\hw4\',...
    'yalefaces_cropped\CroppedYale\yale*']);

Imgs_cr = nan(192 * 168, length(yale_cr_folders));
Imgs_cr_montage = nan(192, 168, length(yale_cr_folders));

% Get one image from the 38 different subjects
for ii = 1:1:length(yale_cr_folders)
    filename = strcat(yale_cr_folders(ii).folder, '\', ...
        yale_cr_folders(ii).name, '\', yale_cr_folders(ii).name, ...
        '_P00A+000E+00.pgm');
    temp = imread(filename);
    Imgs_cr_montage(:, :, ii) = temp;
    Imgs_cr(:, ii) = temp(:);
end

figure, montage(uint8(Imgs_cr_montage))

mean_face_cr = mean(Imgs_cr, 2);
X_Imgs_cr = Imgs_cr - mean_face_cr;

%% SVD
[U_cr, S_cr, V_cr] = svd(X_Imgs_cr, 'econ');
sig_cr = diag(S_cr);

U_cr_rs = nan(size(U_cr));
for jj = 1:1:length(yale_cr_folders)
    U_cr_rs(:, jj) = rescale(U_cr(:, jj), 0, 255);
end
U_cr_rs_mon = reshape(U_cr_rs, [192, 168, length(yale_cr_folders)]);
U_cr_rs_mon = cat(3, U_cr_rs_mon);
figure, montage(uint8(U_cr_rs_mon))

Y_cr = V_cr * X_Imgs_cr';

fig = figure;
fig.Units = 'inches';
fig.Position = [-.1 1.8 6.75 5.0625];
fig.PaperUnits = 'inches';
fig.PaperSize = [6.75 5.0625];

s1 = subplot(2, 2, 1);
s1.Box = 'on';
hold on
h1 = plot(sig_cr, 'ko-', 'LineWidth', 1.1);
h1.MarkerSize = 3.5;
xlabel({'Mode, k'; '(a)'})
ylabel('Singular value, \sigma_k')
axis tight

s2 = subplot(2, 2, 2);
s2.Box = 'on';
hold on
h2 = plot(cumsum(sig_cr) / sum(sig_cr), 'ks-', 'LineWidth', 1.1);
h2.MarkerSize = 3.5;
xlabel({'Mode, k'; '(b)'})
ylabel({'Cumulative energy in', 'first k modes, \sigma_k'})
axis tight

%% LOAD UNCROPPED DATA

addpath('D:\course_work\amath582\hw4\yalefaces_uncropped\yalefaces\')
yale_uncr_folders = dir(['D:\course_work\amath582\hw4\',...
    'yalefaces_uncropped\yalefaces\sub*']);

% Randomly select 38 images (same as cropped dataset) from the dataset
Imgs_uncr = nan(243 * 320, length(yale_cr_folders));
Imgs_uncr_montage = nan(243, 320, length(yale_cr_folders));
rand_files = randperm(length(yale_uncr_folders), length(yale_cr_folders));

c1 = 0;
for ii = rand_files
    c1 = c1 + 1;
    filename = strcat(yale_uncr_folders(ii).folder, '\', ...
        yale_uncr_folders(ii).name);
    temp = imread(filename);
    Imgs_uncr_montage(:, :, c1) = temp;
    Imgs_uncr(:, c1) = temp(:);
end

figure, montage(uint8(Imgs_uncr_montage))

mean_face_uncr = mean(Imgs_uncr, 2);
X_Imgs_uncr = Imgs_uncr - mean_face_uncr;

%% SVD
[U_uncr, S_uncr, V_uncr] = svd(X_Imgs_uncr, 'econ');
sig_uncr = diag(S_uncr);


U_uncr_rs = nan(size(U_uncr));
for jj = 1:1:length(yale_cr_folders)
    U_uncr_rs(:, jj) = rescale(U_uncr(:, jj), 0, 255);
end
U_uncr_rs_mon = reshape(U_uncr_rs, [243, 320, length(yale_cr_folders)]);
U_uncr_rs_mon = cat(3, U_uncr_rs_mon);
figure, montage(uint8(U_uncr_rs_mon))

Y_uncr = V_uncr * X_Imgs_uncr';

figure(fig);
subplot(s1);
h3 = plot(sig_uncr, 'ro--', 'LineWidth', 1.1);
h3.MarkerSize = 3.5;
axis tight
legend([h1, h3], {'Cropped', 'Uncropped'})

subplot(s2);
h4 = plot(cumsum(sig_uncr) / sum(sig_uncr), 'rs--', 'LineWidth', 1.1);
h4.MarkerSize = 3.5;
axis tight
legend([h2, h4], {'Cropped', 'Uncropped'})

%% VISUALIZING THE IMAGES AND SVD MODES

% Mode
cr_imgs = cat(3, uint8(Imgs_cr_montage(:, :, 1:4)),...
    uint8(U_cr_rs_mon(:, :, 1:4)));
uncr_imgs = cat(3, uint8(Imgs_uncr_montage(:, :, 1:4)),...
    uint8(U_uncr_rs_mon(:, :, 1:4)));
fig = figure;
fig.Units = 'inches';
fig.Position = [-.1 1.8 6.75 5.0625];
fig.PaperUnits = 'inches';
fig.PaperSize = [6.75 5.0625];
sa = subplot(2, 1, 1);
montage(cr_imgs, 'Size', [2, 4])
xlabel('(a)')
sb = subplot(2, 1, 2);
montage(uncr_imgs, 'Size', [2, 4])
xlabel('(b)')

% Mode strengths
fig = figure;
fig.Units = 'inches';
fig.Position = [-.1 1.8 6.75 5.0625];
fig.PaperUnits = 'inches';
fig.PaperSize = [6.75 5.0625];
plot_labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)'};
for mm = 1: 1: 3
    sc = subplot(3, 2, (mm - 1) * 2 + 1);
    h = plot(V_cr(:, mm), 'k.-', 'LineWidth', 1.1);
    xlabel(['Mode ', num2str(mm), '      ', plot_labels{(mm - 1) * 2 + 1}])
    ylabel('Magnitude')
    h.MarkerSize = 10;
    sd = subplot(3, 2, (mm - 1) * 2 + 2);
    h = plot(V_uncr(:, mm), 'k.-', 'LineWidth', 1.1);
    xlabel(['Mode ', num2str(mm), '      ', plot_labels{(mm - 1) * 2 + 2}])
    ylabel('Magnitude')
    h.MarkerSize = 10;
end

%% IMAGE RECONSTRUCTION

fig = gobjects([2, 1]);
fig(1) = figure;
fig(2) = figure;
fig(1).Units = 'inches';
fig(1).Position = [-.1 1.8 6.75 6];
fig(1).PaperUnits = 'inches';
fig(1).PaperSize = [6.75 6];
fig(2).Units = 'inches';
fig(2).Position = [-.1 1.8 6.75 6];
fig(2).PaperUnits = 'inches';
fig(2).PaperSize = [6.75 6];
r = [5, 12, 20, 24, 28, 32];
I_recon_cr = nan(192, 168, 32);
I_recon_uncr = nan(243, 320, 32);
for ii = 1: 1: length(r)
    temp_cr = mean_face_cr +...
        U_cr(:, 1:r(ii)) * S_cr(1:r(ii), 1:r(ii)) * V_cr(:, 1:r(ii))';
    temp_uncr = mean_face_uncr +...
        U_uncr(:, 1:r(ii)) * S_uncr(1:r(ii), 1:r(ii)) * V_uncr(:, 1:r(ii))';
    for jj = 1: 1: 4
        I_recon_cr(:, :, (ii - 1) * 4 + jj) = reshape(temp_cr(:, jj), ...
            [192, 168]);
        I_recon_uncr(:, :, (ii - 1) * 4 + jj) = reshape(temp_uncr(:, jj), ...
            [243, 320]);
        figure(fig(1));
        s5 = subplot(3, 2, ii);
        montage(uint8(I_recon_cr(:, :, (ii - 1) * 4 + (1:4))), 'Size', [1, 4])
        xlabel([plot_labels{ii}, '      ', 'Rank = ', num2str(r(ii))])
        figure(fig(2));
        s6 = subplot(3, 2, ii);
        montage(uint8(I_recon_uncr(:, :, (ii - 1) * 4 + (1:2))), 'Size', [1, 2])
        xlabel([plot_labels{ii}, '      ', 'Rank = ', num2str(r(ii))])
    end    
end

%% COMPARING L2 AND FROBENIUS NORMS

l2_cr = nan(37, 1);
l2_uncr = nan(37, 1);
fro_cr = nan(37, 1);
fro_uncr = nan(37, 1);

for ii = 1:1:37
fro_cr(ii) = norm(sig_cr(ii + 1:end));
fro_uncr(ii) = norm(sig_uncr(ii + 1:end));
end
l2_cr_norm = sig_cr(2:end) / sig_cr(1);
l2_uncr_norm = sig_uncr(2:end) / sig_uncr(1);
fro_cr_norm = fro_cr / norm(sig_cr);
fro_uncr_norm = fro_uncr / norm(sig_uncr);

figure(fig)
s3 = subplot(2, 2, 3);
s3.Box = 'on';
h5 = plot(l2_cr_norm, 'k', 'LineWidth', 1.1);
hold on
h6 = plot(l2_uncr_norm, 'r--', 'LineWidth', 1.1);
xlabel({'Mode, k'; '(c)'})
ylabel('Normalized l_2 norm, ||X - X_k|| / ||X||')
axis tight
legend([h5, h6], {'Cropped', 'Uncropped'})


s4 = subplot(2, 2, 4);
s4.Box = 'on';
hold on
h7 = plot(fro_cr_norm, 'k', 'LineWidth', 1.1);
h8 = plot(fro_uncr_norm, 'r--', 'LineWidth', 1.1);
xlabel({'Mode, k'; '(d)'})
ylabel({'Normalized Frobenius norm,', '||X - X_k||_F / ||X||_F'})
axis tight
legend([h7, h8], {'Cropped', 'Uncropped'})
