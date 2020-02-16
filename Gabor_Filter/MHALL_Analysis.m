
clear; close all; clc

%% LOADING AND DEFINING DATASETS

tr_piano = 16;  % recording time in seconds
y_piano = audioread('music1.wav');
L_piano = length(y_piano);
t_piano = (1:L_piano)';
Fs_piano = L_piano / tr_piano;
% p8 = audioplayer(y_piano, Fs_piano); playblocking(p8);

tr_rec = 14;  % recording time in seconds
y_rec = audioread('music2.wav');
L_rec = length(y_rec);
t_rec = (1:L_rec)';
Fs_rec = L_rec / tr_rec;
% p8 = audioplayer(y_rec, Fs_rec); playblocking(p8);

%% FFTs

freq_piano = Fs_piano * ((-L_piano / 2):1:(L_piano / 2 - 1)) / L_piano;
freq_rec = Fs_rec * ((-L_rec / 2):1:(L_rec / 2 - 1)) / L_rec;
y_pianot = fft(y_piano);
y_rect = fft(y_rec);
Y_piano = abs(y_pianot) * 2 / L_piano;
Y_rec = abs(y_rect) * 2 / L_rec;
Y_piano(1) = Y_piano(1) / 2;
Y_rec(1) = Y_rec(1) / 2;
% psd_piano = (1 / (Fs_piano * L_piano)) * abs(Y_piano) .^ 2;
% psd_rec = (1 / (Fs_rec * L_rec)) * abs(Y_rec) .^ 2;

fig1 = figure;
fig1.Units = 'inches';
fig1.Position = [-.1 1.8 6.75 5.0625];
fig1.PaperUnits = 'inches';
fig1.PaperSize = [6.75 5.0625];
s1 = subplot(2,1,1); % Time domain
plot(t_piano / Fs_piano, y_piano, 'k')
hold on
plot(t_rec / Fs_rec, y_rec, 'r-.');
set(gca,'Fontsize', 12), xlabel({'Time, t [s]', '(a)'}), ylabel('y(t)')
axis('tight')
legend('Piano', 'Recorder')

s2 = subplot(2,1,2); % Fourier domain 
plot(freq_piano, fftshift(Y_piano) / max(Y_piano), 'k');
hold on
plot(freq_rec, fftshift(Y_rec) / max(Y_rec), 'r-.');
set(gca, 'Fontsize', 12)
xlabel({'Frequency, \nu [Hz]', '(b)'}), ylabel('Normalized FFT(y)')
ylim([0 1.1])
% xlim([0, max([Fs_piano / 2, Fs_rec / 2])])
xlim([0 5e3])
legend('Piano', 'Recorder')

% print('mhall', '-depsc', '-r600')
% print('mhall', '-dpng', '-r600')

%% SPECTROGRAM WITH GABOR WINDOW AND FILTERING

tslide_piano = 0:1000:L_piano;
tslide_rec = 0:5e3:L_rec;

% Shannon windows
filt_piano = zeros(length(Y_piano), 1);
filt_piano(abs(freq_piano) > 230 & abs(freq_piano) < 350) = 1;
filt_rec = zeros(length(Y_rec), 1);
filt_rec(abs(freq_rec) > 785 & abs(freq_rec) < 1065) = 1;

y_pianoft = y_pianot .* ifftshift(filt_piano);
y_pianof = ifft(y_pianoft);
y_recft = y_rect .* ifftshift(filt_rec);
y_recf = ifft(y_recft);

% Unfiltered Piano Spectrogram
ygt_spec_piano = nan(length(tslide_piano), L_piano);
for ii = 1:length(tslide_piano)
    g_filt_piano = exp(-(t_piano - tslide_piano(ii)) .^ 2 / (2 * 1500 ^ 2)); % Gabor 
    yg_piano = g_filt_piano .* y_piano;
    ygt_piano = fft(yg_piano);
    ygt_spec_piano(ii, :) = abs(fftshift(ygt_piano'));
end

fig = figure;
fig.Units = 'inches';
fig.Position = [-.1 1.8 6.75 5.0625];
fig.PaperUnits = 'inches';
fig.PaperSize = [6.75 5.0625];
s1 = subplot(2, 2, 1);
pcolor(tslide_piano / Fs_piano, freq_piano, ygt_spec_piano.'), 
shading interp
% ylim([0 2500])
xlabel({'Time, t [s]', '(a)'})
ylabel('Frequency, \nu [Hz]')
ylim([0 4000])
s1.FontSize = 12;
colormap(pink)
clearvars ygt_spec_piano

% Filtered Piano spectrogram
ygt_spec_pianof = nan(length(tslide_piano), L_piano);
for ii = 1:length(tslide_piano)
    g_filt_piano = exp(-(t_piano - tslide_piano(ii)) .^ 2 / (2 * 1500 ^ 2)); % Gabor  
    yg_pianof = g_filt_piano .* y_pianof;
    ygt_pianof = fft(yg_pianof);
    ygt_spec_pianof(ii, :) = abs(fftshift(ygt_pianof'));
end

s3 = subplot(2, 2, 3);
pcolor(tslide_piano / Fs_piano, freq_piano, ygt_spec_pianof.'), 
shading interp
% ylim([0 2500])
xlabel({'Time, t [s]', '(c)'})
ylabel('Frequency, \nu [Hz]')
ylim([200 400])
s3.FontSize = 12;
colormap(pink)
clearvars ygt_spec_pianof

% Unfiltered Recorder spectrogram
ygt_spec_rec = nan(length(tslide_rec), L_rec);
for ii = 1:length(tslide_rec)
    g_filt_rec = exp(-(t_rec - tslide_rec(ii)) .^ 2 / (2 * 1200 ^ 2)); % Gabor 
    yg_rec = g_filt_rec .* y_rec;  
    ygt_rec = fft(yg_rec); 
    ygt_spec_rec(ii, :) = abs(fftshift(ygt_rec'));
end

s2 = subplot(2, 2, 2);
pcolor(tslide_rec / Fs_rec, freq_rec, ygt_spec_rec.'), 
shading interp
% ylim([0 2500])
xlabel({'Time, t [s]', '(b)'})
ylabel('Frequency, \nu [Hz]')
ylim([0 8e3])
s2.FontSize = 12;
colormap(pink)
clearvars ygt_spec_rec

% Filter Recorder spectrogram

ygt_spec_recf = nan(length(tslide_rec), L_rec);
for ii = 1:length(tslide_rec)
    g_filt_rec = exp(-(t_rec - tslide_rec(ii)) .^ 2 / (2 * 1200 ^ 2)); % Gabor     
    yg_recf = g_filt_rec .* y_recf;  
    ygt_recf = fft(yg_recf); 
    ygt_spec_recf(ii, :) = abs(fftshift(ygt_recf'));
end

s4 = subplot(2, 2, 4);
pcolor(tslide_rec / Fs_rec, freq_rec, ygt_spec_recf.'), 
shading interp
% ylim([0 2500])
xlabel({'Time, t [s]', '(d)'})
ylabel('Frequency, \nu [Hz]')
s4.YLim = [750 1150];
s4.FontSize = 12;
colormap(pink)
% print('mhall_sp_analysis', '-depsc', '-r600')
% print('mhall_sp_analysis', '-dpng', '-r600')

clearvars ygt_spec_recf
