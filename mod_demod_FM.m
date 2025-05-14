% mod_demod_FM.m
% 描述: 使用 Test.m4a 实现窄带频率调制 (NBFM) 及解调。
%       系统采样频率固定为原始音频采样频率 (例如 48000 Hz)。
%       载波频率 (fc) 和峰值频偏 (freq_dev) 已相应调低。
%       利用 MATLAB 内置函数 fmmod 和 fmdemod。
%       音频文件将保存为 .wav 格式。
%       增强了音频播放和保存部分的鲁棒性。
%       确保 Test.m4a 文件与此脚本在同一目录下，或提供完整路径。

clear; clc; close all;

%% 1. 加载音频与预处理
audio_filename = 'Test.m4a'; % 您的音频文件名
if ~isfile(audio_filename)
    error('音频文件 "%s" 未找到。请确保文件路径正确。', audio_filename);
end
[audio_signal, fs_audio] = audioread(audio_filename);

% 转为单声道
if size(audio_signal, 2) > 1
    audio_signal_mono = mean(audio_signal, 2);
    fprintf('输入音频为立体声，已转换为单声道。\n');
else
    audio_signal_mono = audio_signal;
    fprintf('输入音频为单声道。\n');
end

% 归一化到 [-1, 1]
if max(abs(audio_signal_mono)) > 1e-9
    audio_signal_norm = audio_signal_mono / max(abs(audio_signal_mono));
else
    audio_signal_norm = audio_signal_mono;
end
original_audio = audio_signal_norm; % 在此脚本中，original_audio 将直接用作基带信号

fprintf('原始音频采样频率 (fs_audio): %d Hz\n', fs_audio);
fprintf('原始音频时长: %.2f 秒\n', length(original_audio)/fs_audio);

%% 2. 定义调制与系统参数
% --- 修改核心：固定系统采样频率，并调整fc和freq_dev ---
fs_system = fs_audio; % 将系统采样频率固定为原始音频的采样频率
fprintf('系统采样频率 (fs_system) 固定为: %d Hz\n', fs_system);

BW_audio = 7500; % Hz, 估计的音频带宽
fprintf('估计的音频带宽 (BW_audio): %d Hz\n', BW_audio);

% FM 调制参数 (针对固定的 fs_system 进行调整，实现NBFM)
% 需要满足: fc + freq_dev + BW_audio < fs_system / 2
% 例如: fs_system = 48000 Hz, fs_system/2 = 24000 Hz
% fc + freq_dev < 24000 - 7500 = 16500 Hz
fc = 10000;      % <<--- 载波频率 (Hz), 大幅降低
freq_dev = 5000; % <<--- 峰值频偏 (Hz), 大幅降低
fprintf('载波频率 (fc) 设置为: %d Hz (NBFM)\n', fc);
fprintf('峰值频偏 (freq_dev) 设置为: %d Hz (NBFM)\n', freq_dev);

% 检查参数是否满足条件
if (fc + freq_dev + BW_audio >= fs_system / 2)
    warning('警告: 当前 fc, freq_dev, BW_audio 和 fs_system 的选择可能导致混叠！');
end

BW_FM_Carson = 2 * (freq_dev + BW_audio);
fprintf('卡森法则估算NBFM信号带宽 (BW_FM_Carson): %d Hz\n', BW_FM_Carson);

% 由于 fs_system 已固定为 fs_audio，不再需要上采样原始音频
audio_formodulation = original_audio; % 直接使用原始（归一化后）音频

N = length(audio_formodulation);
t_system = (0:N-1)'/fs_system; % 时间轴与原始音频一致

%% 3. 频率调制 (FM)
% fmmod(x, Fc, Fs, freqdev)
% Fs 现在是 fs_system (即 fs_audio)
fm_modulated_signal = fmmod(audio_formodulation, fc, fs_system, freq_dev);

%% 4. 模拟信道 - 添加高斯白噪声 (可选)
add_noise = 0; % true: 添加噪声; false: 不加噪声
if add_noise
    SNR_dB = 25; % 信噪比 (dB) - NBFM的抗噪性可能不如WBFM
    fprintf('向已调NBFM信号添加 SNR = %d dB 的高斯白噪声。\n', SNR_dB);
    fm_modulated_signal_channel = awgn(fm_modulated_signal, SNR_dB, 'measured');
else
    fprintf('信道为理想信道，未添加噪声。\n');
    fm_modulated_signal_channel = fm_modulated_signal;
end

%% 5. 频率解调 (FM)
% fmdemod(y, Fc, Fs, freqdev)
% Fs 现在是 fs_system (即 fs_audio)

% 失真情况
fc = fc * 0.95;

demodulated_signal_raw = fmdemod(fm_modulated_signal_channel, fc, fs_system, freq_dev);

if max(abs(demodulated_signal_raw)) > 1e-9
    demodulated_signal_normalized = demodulated_signal_raw / max(abs(demodulated_signal_raw));
else
    demodulated_signal_normalized = demodulated_signal_raw;
end

% 由于 fs_system == fs_audio，不再需要下采样
demodulated_signal_output = demodulated_signal_normalized;
fprintf('解调后信号的采样率与原始音频一致: %d Hz。\n', fs_audio);


%% 6. 结果对比 - 播放与绘图
if add_noise
    noise_condition_str = '有';
    playback_noise_desc = '有噪声';
else
    noise_condition_str = '无';
    playback_noise_desc = '无噪声';
end

fprintf('准备播放原始音频...\n');
try
    if ~isempty(original_audio) && all(isfinite(original_audio(:))) && fs_audio > 0
        soundsc(original_audio, fs_audio);
        pause(max(0.1, length(original_audio)/fs_audio + 0.5));
    else
        fprintf('警告: 原始音频数据无效，无法播放。\n');
    end
catch ME_sound_orig
    fprintf('警告: 播放原始音频时发生错误: %s\n', ME_sound_orig.message);
    fprintf('请检查您的音频设备和驱动程序。\n');
end

fprintf('准备播放已调NBFM信号 (%s) (注意: 可能包含高频成分)...\n', playback_noise_desc);
try
    if ~isempty(fm_modulated_signal_channel) && all(isfinite(fm_modulated_signal_channel(:))) && fs_system > 0
        soundsc(fm_modulated_signal_channel, fs_system);
        pause(max(0.1, length(fm_modulated_signal_channel)/fs_system + 0.5));
    else
        fprintf('警告: 已调NBFM信号数据无效，无法播放。\n');
    end
catch ME_sound_mod
    fprintf('警告: 播放已调NBFM信号时发生错误: %s\n', ME_sound_mod.message);
    fprintf('请检查您的音频设备和驱动程序。\n');
end

fprintf('准备播放解调后的NBFM音频 (%s)...\n', playback_noise_desc);
try
    if ~isempty(demodulated_signal_output) && all(isfinite(demodulated_signal_output(:))) && fs_audio > 0
        soundsc(demodulated_signal_output, fs_audio);
        pause(max(0.1, length(demodulated_signal_output)/fs_audio + 0.5));
    else
        fprintf('警告: 解调后NBFM信号数据无效，无法播放。\n');
    end
catch ME_sound_demod
    fprintf('警告: 播放解调后NBFM信号时发生错误: %s\n', ME_sound_demod.message);
    fprintf('请检查您的音频设备和驱动程序。\n');
end

% --- 绘图 ---
t_orig_audio = (0:length(original_audio)-1)' / fs_audio; % 时间轴与 t_system 相同
figure('Name', 'NBFM 调制与解调过程 (fs固定)');

subplot(3,2,1);
plot(t_orig_audio, original_audio);
title('原始音频信号 (时域)'); xlabel('时间 (s)'); ylabel('幅度'); grid on; xlim([0 min(5, t_orig_audio(end))]);

[P1_orig, f_orig] = calculate_spectrum(original_audio, fs_audio);
subplot(3,2,2);
plot(f_orig, P1_orig);
title('原始音频信号 (频域)'); xlabel('频率 (Hz)'); ylabel('幅度谱'); grid on; xlim([0 BW_audio*1.5]);

subplot(3,2,3);
plot(t_system, fm_modulated_signal_channel);
title(sprintf('已调NBFM信号 (%s噪声, fc=%dHz)', noise_condition_str, round(fc)));
xlabel('时间 (s)'); ylabel('幅度'); grid on; xlim([0 min(5, t_orig_audio(end))]); % 时间轴与原始一致

[P1_mod, f_mod] = calculate_spectrum(fm_modulated_signal_channel, fs_system);
subplot(3,2,4);
plot(f_mod, P1_mod);
title(sprintf('已调NBFM信号 (频域, %s噪声)', noise_condition_str));
xlabel('频率 (Hz)'); ylabel('幅度谱'); grid on; xlim([max(0, fc - BW_FM_Carson/2 * 1.2) fc + BW_FM_Carson/2 * 1.2]);

% t_demod_audio 与 t_orig_audio 相同
subplot(3,2,5);
plot(t_orig_audio, demodulated_signal_output);
title(sprintf('解调后NBFM信号 (%s噪声)', noise_condition_str));
xlabel('时间 (s)'); ylabel('幅度'); grid on; xlim([0 min(5, t_orig_audio(end))]);

[P1_demod, f_demod] = calculate_spectrum(demodulated_signal_output, fs_audio);
subplot(3,2,6);
plot(f_demod, P1_demod);
title(sprintf('解调后NBFM信号 (频域, %s噪声)', noise_condition_str));
xlabel('频率 (Hz)'); ylabel('幅度谱'); grid on; xlim([0 BW_audio*1.5]);

sgtitle(sprintf('NBFM 调制与解调 (fs=%dHz, fc=%dHz, fdev=%dHz)', fs_system, round(fc), round(freq_dev)));

%% 7. 保存音频文件 (格式为 .wav)
modulated_output_filename = sprintf('NBFM_mod_Test_fc%dk_fdev%dk.wav', round(fc/1000), round(freq_dev/1000));
demodulated_output_filename = sprintf('NBFM_demod_Test_fc%dk_fdev%dk.wav', round(fc/1000), round(freq_dev/1000));

if ~isempty(fm_modulated_signal_channel) && all(isfinite(fm_modulated_signal_channel(:)))
    try
        audiowrite(modulated_output_filename, fm_modulated_signal_channel, fs_system);
        fprintf('已调NBFM信号已保存到: %s (采样率: %d Hz)\n', modulated_output_filename, fs_system);
    catch ME_write_mod
        fprintf('警告: 保存已调NBFM信号时发生错误: %s\n', ME_write_mod.message);
    end
else
    fprintf('警告: 已调NBFM信号为空或包含无效值，未保存。\n');
end

if ~isempty(demodulated_signal_output) && all(isfinite(demodulated_signal_output(:)))
    try
        audiowrite(demodulated_output_filename, demodulated_signal_output, fs_audio);
        fprintf('解调后NBFM信号已保存到: %s (采样率: %d Hz)\n', demodulated_output_filename, fs_audio);
    catch ME_write_demod
        fprintf('警告: 保存解调后NBFM信号时发生错误: %s\n', ME_write_demod.message);
    end
else
    fprintf('警告: 解调后NBFM信号为空或包含无效值，未保存。\n');
end

fprintf('脚本执行完毕。\n');

%% 8.计算MSE
% 假设 original_audio_norm 是原始归一化单声道信号
% demodulated_signal_output 是最终的解调输出信号
% 确保两者长度一致，如果因为滤波器延迟等原因不一致，需要对齐或截取可比较部分
if length(original_audio) == length(demodulated_signal_output)
    mse_value = mean((original_audio - demodulated_signal_output).^2);
    fprintf('均方误差 (MSE): %e\n', mse_value);
else
    fprintf('警告: 原始信号与解调信号长度不一致，无法直接计算MSE。\n');
    % 可以考虑截取较短的长度进行比较
    % len = min(length(original_audio_norm), length(demodulated_signal_output));
    % mse_value = mean((original_audio_norm(1:len) - demodulated_signal_output(1:len)).^2);
    % fprintf('部分信号的均方误差 (MSE): %e\n', mse_value);
end
%% 辅助函数：计算单边幅度谱
function [P1, f] = calculate_spectrum(signal_data, Fs)
    L = length(signal_data);
    if L == 0 || Fs <= 0
        P1 = [];
        f = [];
        return;
    end
    Y = fft(signal_data);
    P2 = abs(Y/L);
    if mod(L,2) == 0
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
    else
        P1 = P2(1:(L+1)/2);
        P1(2:end) = 2*P1(2:end);
    end
    f = Fs*(0:(L/2))/L;
end
