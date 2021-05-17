% soundFiles = ['2.4kHz to 8kHz'; '100Hz to 300Hz'; '300Hz to 2.4kHz'; 'adult-female-average'; 'adult-female-noisy'; 'adult-female-quiet'; 
% 'Electro_Techno';'Multi_male_speaker';'Opera_Orchestra';'Single_male_average';'Single_male_loud';'Single_male_quiet']

% Pass in audio file name without .mp3 extension
process('Single_male_quiet');

function process(name)
    filename = append(name,'.mp3');

% PHASE 1
% 3.1 Audio Read and find Sampling rate
[y,Fs] = audioread(filename);
info = audioinfo(filename);
disp(info.SampleRate);

% 3.2 Check Size. Add columns if stereo.
[n,m] = size(filename);
if n==2
    y = sum(y,2)/size(y,2);
end

% 3.3 Play sound in Matlab
% sound(y,Fs);

% 3.4 Write the sound to a new file
audiowrite(append(name,'.wav'),y,Fs);

% 3.5.Plot sound waveform as function of sample number
% t = (0:1/Fs:info.Duration)';
% t = t(1:end-1);
% 
% figure(1)
%  plot(t,y)
% xlabel('Sample Number')
% ylabel('Audio Signal')

% 3.6 Downsample to 16KHz if needed. If less, redo 3.1??
if Fs > 16000
    resampled_audio = resample(y, 16000, info.SampleRate);
else
    disp('smaller than 16kHz')
    resampled_audio = y;
end

% PHASE 2
% TASK 5. Filter sounds with the passband bank
fs=16000; %% sampling frequency
data = resampled_audio; %here data is the signal that will be passed through the filter bank

% INITAL DESIGN (Phase 2)
% lowFreqRanges = [100, 200; 201, 500];
% midFreqRanges = [501, 750; 751, 950; 951, 1100; 1101, 1350; 1351, 1800; 1801, 2300; 2301, 2650; 2651, 3300; 3301, 4500;];
% highFreqRanges = [4501, 6000; 6001, 7050];

% REDESIGN TESTING (Phase 3)
% Task 15: SPACING CHANGES
    % Original Spacing + 50 Hz overlap
        % lowFreqRanges = [100, 200; 151, 500];
        % midFreqRanges = [451, 750; 701, 950; 901, 1100; 1051, 1350; 1301, 1800; 1751, 2300; 2251, 2650; 2601, 3300; 3251, 4500;];
        % highFreqRanges = [4451, 6000; 5951, 7050];
    % More Dense space from 800-2000 Hz
        % lowFreqRanges = [100, 200; 201, 300;301,400;401,500;];
        % midFreqRanges = [501, 750; 751, 850; 851, 950; 951, 1150; 1151, 1250; 1251, 1450; 1451, 1650; 1651, 1850; 1851, 2050;];
        % highFreqRanges = [2050, 3500; 3501, 7050];
    % Linear Spacing-13 channels
        % lowFreqRanges = [100, 605; 606, 1211];
        % midFreqRanges = [1212, 1817; 1818, 2422; 2423, 3027; 3028, 3632; 3633, 4238; 4239, 4844; 4845, 5450; 5451, 6056; 6057, 6672;];
        % highFreqRanges = [6673, 7278; 7279, 7884];
    %  Linear Spacing + 100 Hz Overlap -13 channels
        % lowFreqRanges = [100, 605; 506, 1211];
        % midFreqRanges = [1112, 1817; 1718, 2422; 2323, 3027; 2928, 3632; 3533, 4238; 4139, 4844; 4745, 5450; 5351, 6056; 5957, 6572;];
        % highFreqRanges = [6573, 7278; 7179, 7884];
    %  Linear Spacing + 200 Hz Overlap -13 channels
        % lowFreqRanges = [100, 605; 406, 1211];
        % midFreqRanges = [1012, 1817; 1618, 2422; 2223, 3027; 2828, 3632; 3433, 4238; 4039, 4844; 4645, 5450; 5251, 6056; 5857, 6572;];
        % highFreqRanges = [6473, 7278; 7079, 7884];
% Task 16: NUMBER OF CHANNELS CHANGES
    % 1 Channel (similar distribution) - much more noise
        % lowFreqRanges = [100, 500];
        % midFreqRanges = [501, 4500;];
        % highFreqRanges = [4501, 7050];
    %  7 Channels (similar distribution) - more noise
        % lowFreqRanges = [100, 500];
        % midFreqRanges = [501, 750; 751, 1100; 1101, 1800; 1801, 2650; 2651, 4500;];
        % highFreqRanges = [4501, 7050];
% FINAL DESIGN (Phase 3)
% 16 Chanels (simmilar distribution) - similar sound
    lowFreqRanges = [100, 230; 230, 360; 360,500;];
    midFreqRanges = [500, 900; 900, 1300; 1300, 1700; 1700, 2100; 2100, 2500; 2500, 2900; 2900, 3300; 3300, 3700; 3700, 4100; 4100, 4500];
    highFreqRanges = [4500, 5665; 5665, 6833; 6833, 7999];
    %  22 Channels (simmilar distribution)- extra noise at end
        % lowFreqRanges = [100, 200; 200, 300; 300, 400; 400, 500];
        % midFreqRanges = [500, 865; 865, 1230; 1230, 1595; 1595, 1960; 1960, 2325; 2325, 2690; 2690, 3055; 3055, 3420; 3420, 3785; 2785, 4150; 4150, 4500];
        % highFreqRanges = [4500, 5375; 5375, 6350; 6350, 7125; 7125, 8000];
    %  26 Channels (same distibution) - A lot of static at beginning
        % lowFreqRanges = [100, 150; 151, 200; 201, 3510; 351, 500];
        % midFreqRanges = [501, 625; 625, 750; 751, 850; 851, 950; 951, 1025; 1026, 1100; 1101, 1225; 1226, 1350; 1351, 1575; 1576, 1800; 1801, 2050; 2050, 2300; 2301, 2475; 2475, 2650; 2651, 2975; 2976, 3300; 3301, 3900; 3901, 4500;];
        % highFreqRanges = [4501, 5250; 5251, 6000; 6001, 6525; 6526, 7050];
    %  (even distribution of 400Hz per channel) - very noisy 
        % lowFreqRanges = [100, 500];
        % midFreqRanges = [500, 900; 900, 1300; 1300, 1700; 1700, 2100; 2100, 2500; 2500, 2900; 2900, 3300; 3300, 3700; 3700, 4100; 4100, 4500; 4500, 5000];
        % highFreqRanges = [4900, 5200; 5200, 5700; 5700, 6100; 6100, 6500; 6500, 6900; 6900, 7300; 7300, 7700];

set_apass = 1;
set_astop = 80;

% BUTTERWORTH LOWPASS FILTER FOR ENVELOPE
% All frequency values are in Hz.
sampleFreq = 8000;  % Sampling Frequency

N  = 6; % Order
Fc = 400;   % Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass('N,F3dB', N, Fc, sampleFreq);
ButterWorthLowPass = design(h, 'butter');


%LOW FREQUENCIES (Phase 2)
lowFrequencyChannelsEnveloped = zeros(2, length(data));
lowFrequencyChannels = zeros(2, length(data));
for i = 1:size(lowFreqRanges,1)
    % ORIGINAL DESIGN (Phase 2)
    N = 6;    % Order
    Fpass1 = lowFreqRanges(i);  % First Passband Frequency
%   Fpass1_normalized = lowFreqRanges(i)/(fs/2);
    Fpass2 = lowFreqRanges(i, 2);  % Second Passband Frequency
%   Fpass2_normalized = lowFreqRanges(i, 2)/(fs/2);
    Apass  = set_apass;    % Passband Ripple (dB)
%     Astop  = set_astop;   % Stopband Attenuation (dB)  
%     h  = fdesign.bandpass('N,Fp1,Fp2,Ast1,Ap,Ast2', N, Fpass1, Fpass2, Astop, Apass, Astop, fs);
%     Hd = design(h, 'ellip');

    % FINAL REDESIGN (Phase 3)
    % Task 15: FILTER TYPE AND ORDER CHANGES - Use Cheby1 instead of Ellip
        N = 8;    % Order
        Apass = 1;    % Passband Ripple (dB)
        % Construct an FDESIGN object and call its CHEBY1 method.
        h  = fdesign.bandpass('N,Fp1,Fp2,Ap', N, Fpass1, Fpass2, Apass, Fs);
        Hd = design(h, 'cheby1');
        
    % Task 7: Rectify Signal using abs()    
    filtered_audio = abs(filter(Hd, data));
    lowFrequencyChannels(i, :) = transpose(filtered_audio(:,1)); 
    % TASK 8: Envelope channels with Lowpass Filter
    lowFrequencyChannelsEnveloped(i, :) = filter(ButterWorthLowPass, transpose(filtered_audio(:,1)));
end

% TASK 6: Output Lowest Channel
figure
    plot(lowFrequencyChannels(1,:));
    title("Lowest Frequency Channel: 100-200");
    xlabel('Sample Number')
    ylabel('Audio Signal(dB)')
% TASK 9: Output Lowest Enveloped Channel
figure
    plot(lowFrequencyChannelsEnveloped(1,:));
    title("Lowest Frequency Channel: 100-200 (Enveloped)");
    xlabel('Sample Number')
    ylabel('Audio Signal(dB)')
    
% MID FREQUENCIES
midFrequencyChannels = zeros(9, length(data));
midFrequencyChannelsEnveloped = zeros(9, length(data));
for i = 1:size(midFreqRanges, 1)
    N   = 8;    % Order
    Fpass1 = midFreqRanges(i);  % First Passband Frequency
%   Fpass1_normalized = lowFreqRanges(i)/(fs/2);
    Fpass2 = midFreqRanges(i, 2);  % Second Passband Frequency
%   Fpass2_normalized = lowFreqRanges(i, 2)/(fs/2);
    
    % ORIGINAL DESIGN (Phase 2)
    % Construct an FDESIGN object and call its BUTTER method.
%     h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fpass1, Fpass2, fs);
%     Hd = design(h, 'butter');

    % FINAL REDESIGN (Phase 3)
    % Task 15: FILTER TYPE AND ORDER CHANGES - Use Cheby1 instead of Butter
        % Construct an FDESIGN object and call its CHEBY1 method.
        h  = fdesign.bandpass('N,Fp1,Fp2,Ap', N, Fpass1, Fpass2, Apass, Fs);
        Hd = design(h, 'cheby1');
        
%   Task 7: Rectify Signal using abs()
    filtered_audio = abs(filter(Hd, data));
    midFrequencyChannels(i, :) = transpose(filtered_audio(:,1));
%   TASK 8: Envelope channels with Lowpass Filter
    midFrequencyChannelsEnveloped(i, :) = filter(ButterWorthLowPass, transpose(filtered_audio(:,1)));
end

% HIGH FREQUENCIES
highFrequencyChannels = zeros(2, length(data));
highFrequencyChannelsEnveloped = zeros(2, length(data));
for i = 1:size(highFreqRanges,1)    
    N = 6;    % Order
    Fpass1 = highFreqRanges(i);  % First Passband Frequency
%   Fpass1_normalized = lowFreqRanges(i)/(fs/2);
    Fpass2 = highFreqRanges(i, 2);  % Second Passband Frequency
%   Fpass2_normalized = lowFreqRanges(i, 2)/(fs/2);
    Apass  = set_apass;    % Passband Ripple (dB)
    Astop  = set_astop;   % Stopband Attenuation (dB)  
        
    % ORIGINAL DESIGN (Phase 2)
%     h  = fdesign.bandpass('N,Fp1,Fp2,Ast1,Ap,Ast2', N, Fpass1, Fpass2, Astop, Apass, Astop, fs);
%     Hd = design(h, 'ellip');
    
    % FINAL REDESIGN (Phase 3)
        N = 8;    % Order
    % Task 15: FILTER TYPE AND ORDER CHANGES - Use Cheby1 instead of Ellip
    % Construct an FDESIGN object and call its CHEBY1 method.
        h  = fdesign.bandpass('N,Fp1,Fp2,Ap', N, Fpass1, Fpass2, Apass, Fs);
        Hd = design(h, 'cheby1');
        
    % Task 7: Rectify Signal using abs()
    filtered_audio = abs(filter(Hd, data));
    highFrequencyChannels(i, :) = transpose(filtered_audio(:,1));
    % TASK 8: Envelope channels with Lowpass Filter
    highFrequencyChannelsEnveloped(i, :) = filter(ButterWorthLowPass, transpose(filtered_audio(:,1)));

end
% TASK 6: Output Higest Channel
figure
    plot(highFrequencyChannels(2,:));
    title("Highest Frequency Channel: 6001-7050");
    xlabel('Sample Number')
    ylabel('Audio Signal(dB)')
% TASK 9: Output Highest Enveloped Channel
figure
    plot(highFrequencyChannelsEnveloped(2,:));
    title("Highest Frequency Channel: 6001-7050 (Enveloped)");
    xlabel('Sample Number')
    ylabel('Audio Signal(dB)')

% TASK 12: Add signals Together
Fpass1 = lowFreqRanges(1);  % First Passband Frequency
Fpass2 = lowFreqRanges(1, 2);  % Second Passband Frequency
lowCentral = (Fpass1 + Fpass2) /2;
t = 0:length(lowFrequencyChannels(1,:));
x = cos(2*pi*lowCentral*t);

%     TASK 11
y = ammod(x,lowCentral,16000);
y = lowFrequencyChannelsEnveloped(1, :);
outputSignal = y;

% TASK 10: 
for i = 2:size(lowFreqRanges,1)
    Fpass1 = lowFreqRanges(i);  % First Passband Frequency
    Fpass2 = lowFreqRanges(i, 2);  % Second Passband Frequency
    lowCentral = (Fpass1 + Fpass2) /2;
    t = 0:length(lowFrequencyChannels(2,:));
    x = cos(2*pi*lowCentral*t);
    %     TASK 11
    y = ammod(x,lowCentral,16000);
    y = lowFrequencyChannelsEnveloped(i, :);
    outputSignal = outputSignal + y;
end

for i = 1:size(midFreqRanges, 1)
    Fpass1 = midFreqRanges(i);  % First Passband Frequency
    Fpass2 = midFreqRanges(i, 2);  % Second Passband Frequency
    midCentral = (Fpass1 + Fpass2) /2;
    t = 0:length(midFrequencyChannels(2,:));
    x = cos(2*pi*midCentral*t);
    %     TASK 11
    y = ammod(x,midCentral,16000);
    y = midFrequencyChannelsEnveloped(i, :);
    outputSignal = outputSignal + y;
end

for i = 1:size(highFreqRanges,1)
    Fpass1 = highFreqRanges(i);  % First Passband Frequency
    Fpass2 = highFreqRanges(i, 2);  % Second Passband Frequency
    highCentral = (Fpass1 + Fpass2) /2;
    t = 0:length(highFrequencyChannels(2,:));
    x = cos(2*pi*highCentral*t);
%     TASK 11
    y = ammod(x,highCentral,16000);
    y = highFrequencyChannelsEnveloped(i, :);
    outputSignal = outputSignal + y;
end

% normalizedOutput = find(max(outputSignal) == outputSignal);
% TASK 13: 
% sound(outputSignal,16000);
figure
plot(outputSignal)

% PHASE 3
% Normalize Signals to avoid Data Clipping if out of -1 to 1 range
normalizedOutput = outputSignal;
maxValue = abs(max(outputSignal));
disp(maxValue);
if maxValue > 1
    normalizedOutput = outputSignal/maxValue;
end
% Write the filtered sound to a new file
audiowrite(append(name,'_filtered.wav'),normalizedOutput, 16000);
end