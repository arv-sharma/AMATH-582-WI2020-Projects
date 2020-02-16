
clear; close all; clc

%% Load data

load handel
v = y' / 2;

L = length(v);
t = (1:L)';
figure, plot(t / Fs, v);
xlabel('Time [sec]')
ylabel('Amplitude')
title('Signal of Interest, v(n)')

freq_unshift = Fs * (0: (L - 1)) / L;
freq_unshift(((L + 1) / 2 + 1):end) =...
    freq_unshift(((L + 1) / 2 + 1):end) - Fs;   % Frequencies for FFT
freq = fftshift(freq_unshift);
yt = fft(y);
Y = abs(yt) * 2 / L;    % Single-sided spectrum
Y(1) = Y(1) / 2;    % DC adjustment

fig1 = figure;
fig1.Units = 'inches';
fig1.Position = [-.1 1.8 6.75 5.0625];
fig1.PaperUnits = 'inches';
fig1.PaperSize = [6.75 5.0625];
s1 = subplot(2,1,1); % Time domain
plot(t / Fs, y, 'k') 
set(gca,'Fontsize', 12), xlabel({'Time, t [s]', '(a)'}), ylabel('y(t)')
axis('tight')

s2 = subplot(2,1,2); % Fourier domain 
plot(freq, fftshift(Y) / max(Y), 'k');
set(gca, 'Fontsize', 12)
xlabel({'Frequency, \nu [Hz]', '(b)'}), ylabel('Normalized FFT(y)')
axis([0 (Fs / 2) 0 1])

% print('handel', '-depsc', '-r600')
% print('handel', '-dpng', '-r600')

%%  GABOR TRANSFORM - ILLUSTRATION

shannon_filt = zeros(length(t), 1);
[g_filt, g_filt_trans, g_filt_width, shannon_filt(4e4 - 650 : 4e4 + 650)] =...
    deal(exp(-(t - 4e4) .^ 2 / (2 * 400 ^ 2)), exp(-(t - 4.3e4) .^ 2 /...
    (2 * 400 ^ 2)), exp(-(t - 4e4) .^ 2 / (2 * 1000 ^ 2)), 1);  % Different kernels
filters = [g_filt, g_filt_trans, g_filt_width, shannon_filt];
y_filt = filters .* y;
y_filt_t = fft(y_filt);

fig2 = figure(2);
fig2.Units = 'inches';
fig2.Position = [0 4 6.75 5.0625];
fig2.PaperUnits = 'inches';
fig2.PaperSize = [6.75 5.0625];
s1 = subplot(3, 1, 1);
h1 = plot(t / Fs, y, 'k');
hold on 
h2 = plot(t / Fs, g_filt, 'r');
h3 = plot(t / Fs, g_filt_trans, 'm-.');
h4 = plot(t / Fs, g_filt_width, 'g');
h5 = plot(t / Fs, shannon_filt, 'b');
ylabel('y(t), f(t)')
xlabel({'Time, t [s]', '(a)'})
xlim([4.3 5.8])
legend([h1], {'Original signal'})

s2 = subplot(3, 1, 2);
s2.Box = 'on';
hold on
h6 = plot(t / Fs, y_filt(:, 1), 'r');
h7 = plot(t / Fs, y_filt(:, 2), 'm-.');
h8 = plot(t / Fs, y_filt(:, 3), 'g');
h9 = plot(t / Fs, y_filt(:, 4), 'b');
ylabel('y(t)f(t)')
xlabel({'Time, t [s]', '(b)'})
xlim([4.3 5.8])
legend('Gaussian', 'Gaussian, shifted', 'Gaussian, wider', 'Shannon')

s3 = subplot(3,1,3);
s3.Box = 'on';
hold on
h10 = plot(freq, abs(fftshift(y_filt_t(:, 1))) / max(abs(y_filt_t(:, 1))), 'r');
h11 = plot(freq, abs(fftshift(y_filt_t(:, 2))) / max(abs(y_filt_t(:, 2))), 'm-.');
h12 = plot(freq, abs(fftshift(y_filt_t(:, 3))) / max(abs(y_filt_t(:, 3))), 'g');
h13 = plot(freq, abs(fftshift(y_filt_t(:, 4))) / max(abs(y_filt_t(:, 4))), 'b');
xlim([0 1200])
ylim([0 1])
ylabel('FFT(yf)/max(FFT(yf))')
xlabel({'Frequency, \nu [Hz]', '(c)'})

s1.FontSize = 12; s2.FontSize = 12; s3.FontSize = 12;

% print('Gabor_illustration', '-depsc', '-r600')
% print('Gabor_illustration', '-dpng', '-r600')

%% SPECTROGRAM

tslide = 0:500:L;
ygt_spec = nan(length(tslide), L);
for ii = 1:length(tslide)
    g_filt = exp(-(t - tslide(ii)) .^ 2 / (2 * 300 ^ 2)); % Gaussian Gabor kernel
    yg = g_filt .* y; 
    ygt = fft(yg); 
    ygt_spec(ii, :) = abs(fftshift(ygt')); 
%     subplot(3,1,1), plot(t, y, 'k', t, g_filt, 'r')
%     subplot(3,1,2), plot(t, yg, 'k')
%     subplot(3,1,3), plot(freq, abs(fftshift(ygt)) / max(abs(ygt))) 
%     axis([0 2500 0 1])
%     drawnow
%     pause(0.1)
end

fig3 = figure(3);
pcolor(tslide / Fs, freq, ygt_spec.'), 
shading interp
ylim([0 2500])
xlabel('Time, t [s]')
ylabel('Frequency, \nu [Hz]')
% set(gca, 'Ylim', [-50 50], 'Fontsize', [14]) 
colormap(pink)
fig3.Units = 'inches';
fig3.Position = [0 4 6.75 5.0625];
fig3.PaperUnits = 'inches';
fig3.PaperSize = [6.75 5.0625];
ax = gca;
ax.FontSize = 12;
% print('Spectrogram', '-depsc', '-r600')
% print('Spectrogram', '-dpng', '-r600')

%% EFFECT OF CHANGING WINDOW WIDTH

fig4 = figure(4);
labels = {'(a)', '(b)', '(c)', '(d)'};
sig = [50, 200, 1000, 5000];    % Gaussian kernel width parameter
for jj = 1: 1: length(sig)
    
    tslide = 0:150:L;
    ygt_spec = nan(length(tslide), L);
    for ii = 1:length(tslide)
        g_filt = exp(-(t - tslide(ii)) .^ 2 / (2 * sig(jj) ^ 2)); 
        yg = g_filt .* y; 
        ygt = fft(yg); 
        ygt_spec(ii, :) = abs(fftshift(ygt')); 
    end
    
    subplot(2, 2, jj)
    pcolor(tslide / Fs, freq, ygt_spec.'), 
    shading interp
    ylim([0 2500])
    xlabel({'Time, t [s]', labels{jj}})
    ylabel('Frequency, \nu [Hz]')
    set(gca, 'Fontsize', 12) 
    colormap(pink)
end
fig4.Units = 'inches';
fig4.Position = [0 4 6.75 5.0625];
fig4.PaperUnits = 'inches';
fig4.PaperSize = [6.75 5.0625];
% print('Spectrogram_width', '-depsc', '-r600')
% print('Spectrogram_width', '-dpng', '-r600')

%% EFFECT OF CHANGING TRANSLATION PARAMETER

fig5 = figure(5);
slide = [20; 500; 2000; 10000]; % Translation parameter
for jj = 1: 1: length(slide)
    
    tslide = 0:slide(jj):L;
    ygt_spec = nan(length(tslide), L);
    for ii = 1:length(tslide)
        g_filt = exp(-(t - tslide(ii)) .^ 2 / (2 * 200 ^ 2)); 
        yg = g_filt .* y; 
        ygt = fft(yg); 
        ygt_spec(ii, :) = abs(fftshift(ygt')); 
    end
    
    subplot(2, 2, jj)
    pcolor(tslide / Fs, freq, ygt_spec.'), 
    shading interp
    ylim([0 2500])
    xlabel({'Time, t [s]', labels{jj}})
    ylabel('Frequency, \nu [Hz]')
    set(gca, 'Fontsize', 12) 
    colormap(pink) 
end

fig5.Units = 'inches';
fig5.Position = [0 4 6.75 5.0625];
fig5.PaperUnits = 'inches';
fig5.PaperSize = [6.75 5.0625];
% print('Spectrogram_shifted', '-depsc', '-r600')
% print('Spectrogram_shifted', '-dpng', '-r600')

%% EFFECT OF DIFFERENT FILTERS

tslide = 0:500:L;
ygt_spec = nan(length(tslide), L);
ymht_spec = nan(length(tslide), L);
ysht_spec = nan(length(tslide), L);
ytrt_spec = nan(length(tslide), L);

for ii = 1:length(tslide)

    g_filt = exp(-(t - tslide(ii)) .^ 2 / (2 * 200 ^ 2)); % Gaussian filter
    mh_filt = 2 / (sqrt(3 * 200) * pi ^ .25) * ...
            (1 - ((t - tslide(ii)) ./ 200) .^ 2) ...
            .* exp(- (t - tslide(ii)) .^ 2 ./ (2 * 200) ^ 2);
    mh_filt = mh_filt / max(mh_filt);   % Mexican hat filter
    shannon_filt = zeros(length(t), 1); % Shannon filter
    triangle_filt = zeros(length(t), 1);    % Triangular filter
    if tslide(ii) < 650
        shannon_filt(1:tslide(ii) + 650) = 1;
        triangle_filt(1:tslide(ii) + 650) = ...
            interp1([tslide(ii) - 650, tslide(ii), tslide(ii) + 650],...
            [0, 1, 0], 1:tslide(ii) + 650);
    elseif tslide(ii) > (L - 650)
        shannon_filt(tslide(ii) - 650:L) = 1;
        triangle_filt(tslide(ii) - 650:L) = ...
            interp1([tslide(ii) - 650, tslide(ii), tslide(ii) + 650],...
            [0, 1, 0], tslide(ii) - 650:L);    
    else
        shannon_filt(tslide(ii) - 650 : tslide(ii) + 650) = 1;
        triangle_filt(tslide(ii) - 650 : tslide(ii) + 650) = ...
            interp1([tslide(ii) - 650, tslide(ii), tslide(ii) + 650],...
            [0, 1, 0], t(tslide(ii) - 650 : tslide(ii) + 650));
    end
    
    yg = g_filt .* y; 
    ygt = fft(yg); 
    ygt_spec(ii, :) = abs(fftshift(ygt'));
    
    ymh = mh_filt .* y; 
    ymht = fft(ymh); 
    ymht_spec(ii, :) = abs(fftshift(ymht')); 
    
    ysh = shannon_filt .* y; 
    ysht = fft(ysh); 
    ysht_spec(ii, :) = abs(fftshift(ysht'));
    
    ytr = triangle_filt .* y; 
    ytrt = fft(ytr); 
    ytrt_spec(ii, :) = abs(fftshift(ytrt'));
end

fig6 = figure(6);
s1 = subplot(2, 2, 1);
pcolor(tslide / Fs, freq, ygt_spec.'), 
shading interp
ylim([0 2500])
xlabel({'Time, t [s]', '(a)'})
ylabel('Frequency, \nu [Hz]')
colormap(pink)
s2 = subplot(2, 2, 2);
pcolor(tslide / Fs, freq, ymht_spec.'), 
shading interp
ylim([0 2500])
xlabel({'Time, t [s]', '(b)'})
ylabel('Frequency, \nu [Hz]')
colormap(pink)
s3 = subplot(2, 2, 3);
pcolor(tslide / Fs, freq, ysht_spec.'), 
shading interp
ylim([0 2500])
xlabel({'Time, t [s]', '(c)'})
ylabel('Frequency, \nu [Hz]')
colormap(pink)
s4 = subplot(2, 2, 4);
pcolor(tslide / Fs, freq, ytrt_spec.'), 
shading interp
ylim([0 2500])
xlabel({'Time, t [s]', '(d)'})
ylabel('Frequency, \nu [Hz]')
colormap(pink)
s1.FontSize = 12; s2.FontSize = 12; s3.FontSize = 12; s4.FontSize = 12;
fig6.Units = 'inches';
fig6.Position = [0 4 6.75 5.0625];
fig6.PaperUnits = 'inches';
fig6.PaperSize = [6.75 5.0625];
% print('Spectrogram_filters', '-depsc', '-r600')
% print('Spectrogram_filters', '-dpng', '-r600')
