%% LOADING DATASET AND DEFINING DOMAINS

clear; close all; clc;
load Testdata

L = 15; % spatial domain
n = 64; % Fourier modes
x2 = linspace(-L, L, n+1); x = x2(1:n); y = x; z = x;
k = (2 * pi / (2 * L)) * [0:(n / 2 - 1) (-n / 2):-1]; ks = fftshift(k);

[X, Y, Z] = meshgrid(x, y, z);
[Kx, Ky, Kz] = meshgrid(ks, ks, ks);
Unt_composite = nan(n, n, n, 20);

%% AVERAGING THE FF TRANSFORM TO EXTRACT FREQUENCIES OF INTEREST

for j = 1:20
    Un = reshape(Undata(j,:),n,n,n);
    Unt_composite(:, :, :, j) = fftn(Un); % 4-D matrix containing FFT of the 20 signals
end

Unt_ave = mean(Unt_composite, 4); % Average of 20 FFTs
Unts_ave = fftshift(Unt_ave);

[~, maxid] = max(abs(Unts_ave(:))); % Point where magnitude of FFT is maximum

k_x0 = Kx(maxid); k_y0 = Ky(maxid); k_z0 = Kz(maxid);

% Visualization of average FFT in frequency space
close all, isosurface(Kx, Ky, Kz, abs(Unts_ave) / max(abs(Unts_ave(:))), 0.65);
axis([-7 7 -7 7 -7 7]), drawnow
hold on
h2 = plot3(k_x0, k_y0, k_z0, 'r*', 'MarkerSize', 15);
ax = gca; ax.Box = 'on';
xlabel('Wavenumber in x-direction, k_x')
ylabel('Wavenumber in y-direction, k_y') 
zlabel('Wavenumber in z-direction, k_z')
grid on
fig1 = gcf;
fig1.Units = 'inches';
fig1.Position = [-.1 1.8 6.75 5.0625];
fig1.PaperUnits = 'inches';
fig1.PaperSize = [6.75 5.0625];
% print('frequencies', '-depsc', '-r600')
% print('frequencies', '-dpng', '-r600')

fprintf('\nFrequencies of interest are:\nfx = %.3f\nfy = %.3f\nfz = %.3f\n', k_x0, k_y0, k_z0)

%% SIGNAL ANALYSIS

sig = [3, 1, (1 / sqrt(2)), .5]; % Vector with the standard deviations to be used in Gaussian filter
x_int = nan(1, length(sig)); y_int = nan(1, length(sig)); z_int = nan(1, length(sig));

for ii = 1:1:length(sig)
    
    %% FILTER CONSTRUCTION AND VISUALIZATION
    
    % 3-D Gaussian filter
    filter_3 = exp(-((((ifftshift(Kx) - k_x0) .^ 2 ...
        + (ifftshift(Ky) - k_y0) .^ 2 + (ifftshift(Kz) - k_z0) .^ 2)) / (2 * sig(ii) ^ 2)));
    
    % Filter visualization
    close all, isosurface(Kx, Ky, Kz, fftshift(filter_3), .5);
    axis([-7 7 -7 7 -7 7]), drawnow
    xlabel('Kx'), ylabel('Ky'), zlabel('Kz')
    fig2 = gcf;
    pause(1), close

    %% MARBLE PATH ESTIMATION

    x_path = nan(20, 1); y_path = nan(20, 1); z_path = nan(20, 1);
    for jj = 1: 1: 20
        Untf = Unt_composite(:, :, :, jj) .* filter_3;
        Unf = ifftn(Untf); % Inverse FFT of filtered signal in frequency domain

        [~, maxid1] = max(abs(Unf(:)));
        x_path(jj) = X(maxid1); y_path(jj) = Y(maxid1); z_path(jj) = Z(maxid1); % Estimated marble location

    %     close all, isosurface(X, Y, Z, abs(Unf) / max(abs(Unf(:))), 0.9)
    %     axis([-20 20 -20 20 -20 20]), grid on, drawnow
    %     xlabel('x'), ylabel('y'), zlabel('z')
    %     hold on
    %     plot3(x_int(jj), y_int(jj), z_int(jj), 'r*', 'MarkerSize', 20)
    %     pause(1)
    %     close

    end
    
    x_int(ii) = x_path(end); y_int(ii) = y_path(end); z_int(ii) = z_path(end); % Marble location at the final time point

    fprintf('\nIteration %d :', ii)
    fprintf('\nThe marble''s estimated location is:\n')
    fprintf('x = %.3f\ny = %.3f\nz = %.3f\n', x_path(end), y_path(end), z_path(end))

end

% Check to see if the marble location estimates are consistent

if range(x_int) < 1e-5 && range(y_int) < 1e-5 && range(z_int) < 1e-5
    fprintf('\nThe marble''s location where the acoustic beam is to be focused is:\n')
    fprintf('x = %.3f\ny = %.3f\nz = %.3f\n', x_int(end), y_int(end), z_int(end))
else
    warning('The marble location is uncertain. Verify before proceeding.')
end

%% VISUALIZATION OF THE PATH OF THE MARBLE


fig3 = figure;
fig3.Units = 'inches';
fig3.Position = [-.1 1.8 6.75 5.0625];
fig3.PaperUnits = 'inches';
fig3.PaperSize = [6.75 5.0625];
h3 = plot3(x_path, y_path, z_path);
h3.Color = 'k';
h3.LineStyle = ':';
h3.Marker = '.'; h3.MarkerSize = 18;
hold on
h4 = quiver3(x_path(1:19), y_path(1:19), z_path(1:19), diff(x_path), diff(y_path), diff(z_path));
h4.Color = [0.4 0.4 0.4]; h4.MaxHeadSize = 0.2;
h5 = plot3(x_path(end), y_path(end), z_path(end), 'k.', 'MarkerSize', 20);
h5.Marker = 'p'; h5.MarkerFaceColor = 'k'; h5.MarkerSize = 15;
xlabel('x'), ylabel('y'), zlabel('z')
axis([-11 11 -11 11 -11 11]), grid on
ax = gca; ax.Box = 'on'; ax.View = [-35.1 36.8];
% print('path', '-depsc', '-r600')
% print('path', '-dpng', '-r600')
