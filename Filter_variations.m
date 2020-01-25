clear; close all; clc

%% NOISE ON SIGNAL

L = 30; % time slot to transform 
n = 512; % number of Fourier modes 2^9
t2 = linspace(-L, L, n+1); t = t2(1:n); % time discretization 
k = (2 * pi / (2 * L)) * [0:(n / 2 - 1) (-n / 2):-1]; % frequency components of FFT
u = sech(t); % ideal signal in the time domain 

%%  SIGNAL AND SPECTRUM

noise = 10;
ut = fft(u); 
unt = ut + noise * (randn(1, n) + 1i * randn(1, n)); 
un = ifft(unt);

%% GAUSSIAN FILTERS WITH DIFFERENT VARIANCES

colorvec = lines(4);
lnstyles = {'-.', '--', ':', '-'}; 
fig = figure;
fig.Units = 'inches';
fig.Position = [-.1 1.8 6.75 4];
fig.PaperUnits = 'inches';
fig.PaperSize = [6.75 4];
s1 = subplot(3, 1, 1);
s1.Box = 'on';
plot(fftshift(k), abs(fftshift(unt)) / max(abs(fftshift(unt))), 'k', 'LineWidth', 0.5)
hold on
axis([-18 18 0 1])
xlabel({'wavenumber, k', '(a)'}), ylabel('|ut| / max(|ut|)')
s2 = subplot(3, 1, 2);
s2.Box = 'on';
hold on
axis([-18 18 0 1])
xlabel({'wavenumber, k', '(b)'}), ylabel('|ut| / max(|ut|)')
s3 = subplot(3,1,3);
s3.Box = 'on';
plot(t, u, 'k', 'LineWidth', 0.5)
hold on
axis([-18 18 0 1.2])
xlabel({'time, t', '(c)'}), ylabel('|u|')
c1 = 0;
for sig = [3, 1, (1 / sqrt(2)), .5]
    c1 = c1 + 1;
    subplot(s1)
    filter = exp(-(k).^2 / (2 * sig ^ 2)); 
    unft = filter .* unt; 
    unf = ifft(unft);
    plot(fftshift(k),fftshift(filter), 'Color', colorvec(c1, :), ...
        'LineStyle', lnstyles{c1}, 'Linewidth', 1.5)
    
    subplot(s2)
    plot(fftshift(k), abs(fftshift(unft)) / max(abs(fftshift(unft))),...
         'LineStyle', lnstyles{c1}, 'Color', colorvec(c1, :), 'LineWidth', 1.5)

    subplot(s3)
    plot(t, unf, 'Color', colorvec(c1, :), 'LineStyle', lnstyles{c1}, 'Linewidth', 1.5)
end
subplot(s1)
legend('Original signal', '\sigma = 3', '\sigma = 1', '\sigma = \surd 2', '\sigma = 0.5')
subplot(s2)
legend('\sigma = 3', '\sigma = 1', '\sigma = \surd 2', '\sigma = 0.5')
subplot(s3)
legend('Original signal', '\sigma = 3', '\sigma = 1', '\sigma = \surd 2', '\sigma = 0.5')
% print('filter_sig_variations', '-depsc', '-r600')
% print('filter_sig_variations', '-dpng', '-r600')