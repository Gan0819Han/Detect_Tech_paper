% mod_demod_AM.m
% 描述: 使用 Test.m4a 实现双边带大载波幅度调制 (AM DSB-LC)
%       载波频率 fc 设置为 16000 Hz。
%       及包络检波解调，并对比有无信道噪声的情况。
%       增加播放已调信号和保存已调、解调后音频文件的功能。
%       增强了音频播放部分的鲁棒性，使用try-catch处理音频设备错误。
%       确保 Test.m4a 文件与此脚本在同一目录下，或提供完整路径。

clear; clc; close all;

%% 1. 加载音频与预处理
audio_filename = 'Test.m4a'; % 您的音频文件名
if ~isfile(audio_filename)
    error('音频文件 "%s" 未找到。请确保文件路径正确。', audio_filename);
end
[audio_signal, fs_audio] = audioread(audio_filename);

% 转为单声道 (如果原始是立体声)
if size(audio_signal, 2) > 1
    audio_signal_mono = mean(audio_signal, 2);
    fprintf('输入音频为立体声，已转换为单声道。\n');
else
    audio_signal_mono = audio_signal;
    fprintf('输入音频为单声道。\n');
end

% 归一化到 [-1, 1]
if max(abs(audio_signal_mono)) > 1e-9 % 避免除以0或极小值
    audio_signal_norm = audio_signal_mono / max(abs(audio_signal_mono));
else
    audio_signal_norm = audio_signal_mono; % 如果信号几乎为0，则保持原样
end
original_audio = audio_signal_norm; % 保留一份用于对比

fprintf('原始音频采样频率: %d Hz\n', fs_audio);
fprintf('原始音频时长: %.2f 秒\n', length(original_audio)/fs_audio);

%% 2. 定义调制与系统参数
BW_audio = 7500; % Hz, 音频信号带宽的估计值
fprintf('估计的音频带宽 (BW_audio): %d Hz\n', BW_audio);

Ac = 1;            % 载波幅度
fc = 16000;        % <<--- 修改点: 载波频率 (Hz) 设置为 16000 Hz
m = 1.4;           % 调制指数 (m <= 1 避免包络检波失真)
fprintf('载波频率 (fc): %d Hz\n', fc);
fprintf('调制指数 (m): %.2f\n', m);

required_fs_system_nyquist = 2 * (fc + BW_audio);
suggested_fs_system_carrier = 6 * fc; % 经验值，确保载波被良好表示
fs_system = max(required_fs_system_nyquist, suggested_fs_system_carrier);

% 如果原始 fs_audio 已经很高，或者计算出的 fs_system 不比它高多少，调整
if fs_system <= fs_audio * 1.1 % 如果提升不大，至少提升一个有意义的倍数
    % 确保 fs_system 至少是能表示载波的几倍，并且满足奈奎斯特
    fs_system = max(fs_audio * ceil(suggested_fs_system_carrier / fs_audio), required_fs_system_nyquist * 1.2);
    if fs_system < suggested_fs_system_carrier % 再次确保载波表示
         fs_system = suggested_fs_system_carrier;
    end
end
fs_system = round(fs_system); % 取整
fprintf('计算得到的系统采样频率 (fs_system): %d Hz\n', fs_system);

if fs_system == fs_audio
    audio_upsampled = original_audio;
    fprintf('系统采样频率与原始音频采样频率相同，未进行上采样。\n');
elseif fs_system < fs_audio
    warning('计算得到的系统采样率低于原始音频采样率，可能存在参数配置问题。将使用原始采样率进行后续操作。');
    fs_system = fs_audio;
    audio_upsampled = original_audio;
else
    audio_upsampled = resample(original_audio, fs_system, fs_audio);
    fprintf('原始音频已从 %d Hz 上采样到 %d Hz。\n', fs_audio, fs_system);
end

N = length(audio_upsampled);
t = (0:N-1)' / fs_system;

%% 3. 幅度调制 (AM DSB-LC)
carrier = Ac * cos(2*pi*fc*t);
am_modulated_signal = Ac * (1 + m * audio_upsampled) .* carrier;

%% 4. 模拟信道 - 添加高斯白噪声 (可选)
add_noise = 1; % true: 添加噪声 (失真); false: 不加噪声 (理想)
if add_noise
    SNR_dB = 15; % 信噪比 (dB)
    fprintf('向已调信号添加 SNR = %d dB 的高斯白噪声。\n', SNR_dB);
    am_modulated_signal_channel = awgn(am_modulated_signal, SNR_dB, 'measured');
else
    fprintf('信道为理想信道，未添加噪声。\n');
    am_modulated_signal_channel = am_modulated_signal;
end

%% 5. 幅度解调 - 包络检波器
rectified_signal = abs(am_modulated_signal_channel);

lpf_cutoff = BW_audio * 1.1; % LPF截止频率，略高于音频带宽
filter_order = 6; % 滤波器阶数
[b_lpf, a_lpf] = butter(filter_order, lpf_cutoff/(fs_system/2), 'low');
demodulated_signal_raw = filter(b_lpf, a_lpf, rectified_signal);

% 去除直流分量
demodulated_signal_processed = demodulated_signal_raw - mean(demodulated_signal_raw);

% 幅度归一化
if max(abs(demodulated_signal_processed)) > 1e-9 % 避免除以0
    demodulated_signal_normalized = demodulated_signal_processed / max(abs(demodulated_signal_processed));
else
    demodulated_signal_normalized = demodulated_signal_processed; % 若信号几乎为0
end

% 下采样回原始音频采样率
if fs_system ~= fs_audio
    demodulated_signal_output = resample(demodulated_signal_normalized, fs_audio, fs_system);
    fprintf('解调后信号已从 %d Hz 下采样到 %d Hz。\n', fs_system, fs_audio);
else
    demodulated_signal_output = demodulated_signal_normalized;
end

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
        pause(max(0.1, length(original_audio)/fs_audio + 0.5)); % 确保 pause 时长至少为0.1s
    else
        fprintf('警告: 原始音频数据无效，无法播放。\n');
    end
catch ME_sound_orig
    fprintf('警告: 播放原始音频时发生错误: %s\n', ME_sound_orig.message);
    fprintf('请检查您的音频设备和驱动程序。\n');
end

fprintf('准备播放已调AM信号 (%s) (注意: 包含高频载波)...\n', playback_noise_desc);
try
    if ~isempty(am_modulated_signal_channel) && all(isfinite(am_modulated_signal_channel(:))) && fs_system > 0
        soundsc(am_modulated_signal_channel, fs_system);
        pause(max(0.1, length(am_modulated_signal_channel)/fs_system + 0.5));
    else
        fprintf('警告: 已调AM信号数据无效，无法播放。\n');
    end
catch ME_sound_mod
    fprintf('警告: 播放已调AM信号时发生错误: %s\n', ME_sound_mod.message);
    fprintf('请检查您的音频设备和驱动程序。\n');
end

fprintf('准备播放解调后的AM音频 (%s)...\n', playback_noise_desc);
try
    if ~isempty(demodulated_signal_output) && all(isfinite(demodulated_signal_output(:))) && fs_audio > 0
        soundsc(demodulated_signal_output, fs_audio);
        pause(max(0.1, length(demodulated_signal_output)/fs_audio + 0.5));
    else
        fprintf('警告: 解调后AM信号数据无效，无法播放。\n');
    end
catch ME_sound_demod
    fprintf('警告: 播放解调后AM信号时发生错误: %s\n', ME_sound_demod.message);
    fprintf('请检查您的音频设备和驱动程序。\n');
end

% --- 绘图 ---
t_orig_audio = (0:length(original_audio)-1)' / fs_audio;
figure('Name', 'AM DSB-LC 调制与解调过程 (fc=16000Hz)'); % 更新图表标题

subplot(3,2,1);
plot(t_orig_audio, original_audio);
title('原始音频信号 (时域)'); xlabel('时间 (s)'); ylabel('幅度'); grid on; xlim([0 min(5, t_orig_audio(end))]);

[P1_orig, f_orig] = calculate_spectrum(original_audio, fs_audio);
subplot(3,2,2);
plot(f_orig, P1_orig);
title('原始音频信号 (频域)'); xlabel('频率 (Hz)'); ylabel('幅度谱'); grid on; xlim([0 BW_audio*1.5]);

subplot(3,2,3);
plot(t, am_modulated_signal_channel);
title(sprintf('已调AM信号 (%s噪声, fc=%dHz)', noise_condition_str, round(fc)));
xlabel('时间 (s)'); ylabel('幅度'); grid on; xlim([0 min(5, t_orig_audio(end)) * (fs_audio/fs_system) ]);

[P1_mod, f_mod] = calculate_spectrum(am_modulated_signal_channel, fs_system);
subplot(3,2,4);
plot(f_mod, P1_mod);
title(sprintf('已调AM信号 (频域, %s噪声)', noise_condition_str));
xlabel('频率 (Hz)'); ylabel('幅度谱'); grid on; xlim([max(0, fc-BW_audio*1.5) fc+BW_audio*1.5]);

t_demod_audio = (0:length(demodulated_signal_output)-1)' / fs_audio;
subplot(3,2,5);
plot(t_demod_audio, demodulated_signal_output);
title(sprintf('解调后AM信号 (%s噪声)', noise_condition_str));
xlabel('时间 (s)'); ylabel('幅度'); grid on; xlim([0 min(5, t_demod_audio(end))]);

[P1_demod, f_demod] = calculate_spectrum(demodulated_signal_output, fs_audio);
subplot(3,2,6);
plot(f_demod, P1_demod);
title(sprintf('解调后AM信号 (频域, %s噪声)', noise_condition_str));
xlabel('频率 (Hz)'); ylabel('幅度谱'); grid on; xlim([0 BW_audio*1.5]);

sgtitle(sprintf('AM DSB-LC 调制与包络检波解调 (fc=%dHz)',round(fc))); % 更新图表总标题

%% 7. 保存音频文件
modulated_output_filename = 'AM_mod_Test_fc16k.wav'; % 更新文件名以区分
demodulated_output_filename = 'AM_demod_Test_fc16k.wav'; % 更新文件名以区分

if ~isempty(am_modulated_signal_channel) && all(isfinite(am_modulated_signal_channel(:)))
    try
        audiowrite(modulated_output_filename, am_modulated_signal_channel, fs_system);
        fprintf('已调AM信号已保存到: %s (采样率: %d Hz)\n', modulated_output_filename, fs_system);
    catch ME_write_mod
        fprintf('警告: 保存已调AM信号时发生错误: %s\n', ME_write_mod.message);
    end
else
    fprintf('警告: 已调AM信号为空或包含无效值，未保存。\n');
end

if ~isempty(demodulated_signal_output) && all(isfinite(demodulated_signal_output(:)))
    try
        audiowrite(demodulated_output_filename, demodulated_signal_output, fs_audio);
        fprintf('解调后AM信号已保存到: %s (采样率: %d Hz)\n', demodulated_output_filename, fs_audio);
    catch ME_write_demod
        fprintf('警告: 保存解调后AM信号时发生错误: %s\n', ME_write_demod.message);
    end
else
    fprintf('警告: 解调后AM信号为空或包含无效值，未保存。\n');
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
    if L == 0 || Fs <= 0 % 添加Fs有效性检查
        P1 = [];
        f = [];
        % warning('频谱计算跳过：信号为空或采样率无效。'); % 可以取消注释以显示更多警告
        return;
    end
    Y = fft(signal_data);
    P2 = abs(Y/L);
    if mod(L,2) == 0 % 偶数长度
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
    else % 奇数长度
        P1 = P2(1:(L+1)/2);
        P1(2:end) = 2*P1(2:end);
    end
    f = Fs*(0:(L/2))/L;
end
