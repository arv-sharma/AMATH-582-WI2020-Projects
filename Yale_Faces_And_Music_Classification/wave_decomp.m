%% FUNCTIION TO EXTRACT SAMPLES AND PERFORM WAVELET DECOMPOSITION

function [wavelets, freq] = wave_decomp(filename, num_of_samples, Fs_red)
    wavelets = nan(100 * 5 * Fs_red, num_of_samples);
    freq = nan(100, num_of_samples);
    info = audioinfo(filename);
    Fs = info.SampleRate;
    timing_vec = 1:5:info.Duration - 5;
    sample_times = timing_vec(randperm(length(timing_vec), num_of_samples));
    for kk = 1: num_of_samples
        sample = audioread(filename, [sample_times(kk) * Fs + 1, ...
            sample_times(kk) * Fs + 5 * Fs]);
        sample_mono = mean(sample, 2);
        resampled_mono = resample(sample_mono, Fs_red, Fs);
        [resampled_mono_wav, freq(:, kk)] = cwt(resampled_mono, 'morse',...
            Fs_red, 'WaveletParameters', [750, 750 * 40]);
        wavelets(:, kk) = abs(resampled_mono_wav(:));
    end
    freq = freq(:, 1);
end