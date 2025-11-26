audio_dir = "../audio/";
audio_files = dir(fullfile(audio_dir, "*.wav"));
name = string({audio_files.name});
snr_increment = -0.5;
snr_levels = 5:snr_increment:-5;

% Exclude noise and snr files
is_noise = contains(names, "_noise");
is_snr = contains(names, "_SNR(");

clean_mask = ~is_noise & is_snr;
clean_files = audio_files(clean_mask);

for i = 1:numel(clean_files)
    clean_name = clean_files(i);
    clean_path = fullfile(audio_dir, clean_name);

    % get basename of file
    [~, base, ~] = fileparts(clean_name);

    %noise file : basename_noise.wav
    noise_name = base + "_noise.wav";
    noise_path = fullfile(audio_dir, noise_name);

    if ~isfile(noise_path)
        fprintf("Skipping due to no noise file found:", base, noise_name);
        continue
    end

    fprintf("\nProcessing: %\n", base);

    % Load audio files
    [clean, fs] = audioread(clean_path);
    [noise, fs2] = audioread(noise_path);

    % check sampling rate
    % if fs ~= fs2
    
    % Match lengths to the clean
    L = length(clean);
    if length(noise) >= L
        noise = noise(1:L);
    else
        noise = [noise; zeros(L - length(noise), 1)]; % Pad noise with zeros
    end

    % Power computation
    clean_power = mean(clean.^2);
    noise_power = mean(noise.^2);

    % Now mix signals
    for snr = snr_levels
        % compute target noise to reach SNR
        target_noise_power = clean_power / (10^(snr/10));
        scale = sqrt(target_noise_power/noise_power);

        % Get the noise amplified
        noise_scaled = noise * scale;
        % Mix signals
        mixed_signal = clan + noise_scaled;
        mixed_signal = (max(abs(mixed_signal)))+1e-8);
        
        % get the basename with SNR and save it
        snr_dir = fullfile(audio_dir, "SNRs");
        if ~isfolder(snr_dir)
            mkdir(snr_dir);
        end
        out_name = sprintf("%s_SNR(%+1f).wav", base, snr);
        out_path = (fullfile(snr_dir, out_name));

        % safe file and print what had done
        audiowrite(out_path, mixed_signals, fs);
        fprintf(" Saved: %s\n", out_name);
    end
end
