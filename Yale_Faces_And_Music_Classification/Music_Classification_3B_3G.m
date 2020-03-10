%% MUSIC CLASSIFICATION - 3 BANDS, 3 GENRES

clear; clc; close all

%% BUILDING TRAINING AND TESTING DATASETS

addpath(genpath('D:\course_work\amath582\hw4\music_data\'))

jazz_files = dir(['D:\course_work\amath582\hw4\music_data\',...
    'jazz\keithmitchell_*']);
carnatic_files = dir(['D:\course_work\amath582\hw4\music_data\',...
    'carnatic\semma*']);
classical_files = dir(['D:\course_work\amath582\hw4\music_data\',...
    'classical\Beethoven*']);

all_files = cell(3, 1);
[all_files{1}, all_files{2}, all_files{3}] = deal(...
    extractfield(jazz_files, 'name'), extractfield(carnatic_files,...
    'name'), extractfield(classical_files, 'name'));

Fs_red = 15e3;  % Reduced sampling rate
% size = (wvlt_decomp * num_of_samples * num_of_files/case, num_of_cases)
training_set = nan(75e5, 144);
test_set = nan(75e5, 36);
for ii = 1: 1: 3
    files = all_files{ii};
    files_order = randperm(length(files));  % Randomly pick files to train on
    for jj = 1: 1: 10        
        filename = files{files_order(jj)};
        if jj <= 8
            training_set(:, (ii - 1) * 48 + (jj - 1) * 6 + (1:6)) =...
                wave_decomp(filename, 6, Fs_red);
        else
            test_set(:, (ii - 1) * 12 + (jj - 9) * 6 + (1:6)) =...
                wave_decomp(filename, 6, Fs_red);
        end
    end
end

%% SVD OF TRAINING DATA

mean_training_set = mean(training_set, 2);
[U, S, V] = svd(training_set - mean_training_set, 'econ');

% figure,
% plot(diag(S), 'ko')

% SVD Modes

[~, freq] = wave_decomp(filename, 1, Fs_red);
fig = figure;
fig.Units = 'inches';
fig.Position = [-.1 1.8 6.75 5.0625];
fig.PaperUnits = 'inches';
fig.PaperSize = [6.75 5.0625];
plot_labels = {'(a)', '(b)', '(c)', '(d)'};
for mm = 1: 1: 4
    s = subplot(2, 2, mm);
    pcolor((1: Fs_red * 5) / Fs_red, freq, reshape(abs(U(:, mm)),...
        [100, Fs_red * 5]))
    shading interp
    axis tight
    colormap(pink)
    xlabel({'Time [s]'; plot_labels{mm}})
    ylabel('Frequency [Hz]')
end

% Projection of data onto principal modes

SVD_proj = (U' * training_set)';

fig = figure;
fig.Units = 'inches';
fig.Position = [-.1 1.8 6.75 6];
fig.PaperUnits = 'inches';
fig.PaperSize = [6.75 6];
plot_labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};
for mm = 1:1:6
    s1 = subplot(6, 1, mm);
    s1.Box = 'on';
    hold on
    h1 = histogram((SVD_proj(1:48, mm)), 'EdgeColor', 'r', 'FaceColor', 'r');
    h2 = histogram((SVD_proj(49:96, mm)), 'EdgeColor', 'm', 'FaceColor', 'm');
    h3 = histogram((SVD_proj(97:144, mm)), 'EdgeColor', 'b', 'FaceColor', 'b');
    h2.LineStyle = '--'; h3.LineStyle = ':';
    h1.LineWidth = 1.1; h2.LineWidth = 1.1; h3.LineWidth = 1.1;
    xlabel(plot_labels{mm})
    ylabel('Counts')
    legend('Jazz', 'Carnatic', 'Classical')
end

fig = figure;
fig.Units = 'inches';
fig.Position = [-.1 1.8 6.75 5.0625];
fig.PaperUnits = 'inches';
fig.PaperSize = [6.75 5.0625];
ax = gca;
ax.Box = 'on';
ax.View = [11.412, 63.0438];
hold on
h1 = plot3(SVD_proj(1:48, 1), SVD_proj(1:48, 2), SVD_proj(1:48, 6), 'ro');
h2 = plot3(SVD_proj(49:96, 1), SVD_proj(49:96, 2), SVD_proj(49:96, 6), 'm^');
h3 = plot3(SVD_proj(97:144, 1), SVD_proj(97:144, 2), SVD_proj(97:144, 6), 'bp');
xlabel('PCA 1'), ylabel('PCA 2'), zlabel('PCA 6')
legend('Jazz', 'Carnatic', 'Classical')

% figure;
% for mm = 1:1:6
% subplot(6, 3, (mm - 1) * 3 + 1)
% plot(1:48, (SVD_proj(1:48, mm)), 'ko-')
% subplot(6, 3, (mm - 1) * 3 + 2)
% plot(1:48, (SVD_proj(49:96, mm)), 'ko-')
% subplot(6, 3, (mm - 1) * 3 + 3)
% plot(1:48, (SVD_proj(97:144, mm)), 'ko-')
% end

%% LDA WITH CROSS VALIDATION

labels = [-1 * ones(40, 1);...
          0 * ones(40, 1);...
          1 * ones(40, 1)];
truth = [-1 * ones(8, 1);...
          0 * ones(8, 1);...
          1 * ones(8, 1)];
accuracy_lin = nan(length(truth), 1);
accuracy_quad = nan(length(truth), 1);
for nn = 1: 1: 50
    rand_vec = randperm(48);
    train_mat = [SVD_proj(rand_vec(1:40), [1, 2, 6]);...
                 SVD_proj(rand_vec(1:40) + 48, [1, 2, 6]);...
                 SVD_proj(rand_vec(1:40) + 96, [1, 2, 6])];
    test_mat = [SVD_proj(rand_vec(41:48), [1, 2, 6]);...
                SVD_proj(rand_vec(41:48) + 48, [1, 2, 6]);...
                SVD_proj(rand_vec(41:48) + 96, [1, 2, 6])];
    class_lin = classify(test_mat, train_mat, labels);
    class_quad = classify(test_mat, train_mat, labels, 'quadratic');
    accuracy_lin(nn, 1) = sum(class_lin == truth) / length(truth);
    accuracy_quad(nn, 1) = sum(class_quad == truth) / length(truth);    
end

fprintf('\nTraining data:\n')
fprintf('Quadratic discriminant accuracy\n')
fprintf('Mean = %.2f \t Minimum = %.2f \t Maximum = %.2f\n',...
    mean(accuracy_quad), min(accuracy_quad), max(accuracy_quad))
fprintf('Linear discriminant accuracy\n')
fprintf('Mean = %.2f \t Minimum = %.2f \t Maximum = %.2f\n',...
    mean(accuracy_lin), min(accuracy_lin), max(accuracy_lin))

fig = figure;
fig.Units = 'inches';
fig.Position = [-.1 1.8 6.75 5.0625];
fig.PaperUnits = 'inches';
fig.PaperSize = [6.75 5.0625];
ax = gca;
ax.Box = 'on';
h1 = bar([accuracy_lin, accuracy_quad]);
legend('Linear', 'Quadratic')
xlabel('Trial')
ylabel('Accuracy')

%% LDA ON TEST DATA

test_truth = [-1 * ones(12, 1);...
              0 * ones(12, 1);...
              1 * ones(12, 1)];
test_mat =(U(:, [1, 2, 6])' * test_set)';
labels = [-1 * ones(48, 1);...
          0 * ones(48, 1);...
          1 * ones(48, 1)];
train_mat = [SVD_proj(1:48, [1, 2, 6]);
             SVD_proj(49:96, [1, 2, 6]);
             SVD_proj(97:144, [1, 2, 6])];
class_lin = classify(test_mat, train_mat, labels);
class_quad = classify(test_mat, train_mat, labels, 'quadratic');
accuracy_lin = sum(class_lin == test_truth) / length(test_truth);
accuracy_quad = sum(class_quad == test_truth) / length(test_truth);    

fprintf('\nTest data:\n')
fprintf('Quadratic discriminant\n')
fprintf('Accuracy = %.2f\n', accuracy_quad)
fprintf('Linear discriminant\n')
fprintf('Accuracy = %.2f\n', accuracy_lin)
