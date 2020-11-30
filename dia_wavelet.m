% Read the saved diameter file
dia_name = uigetfile('*.mat');
load(dia_name);
[dia_row, dia_col] = size(dia_mat);

% Low pass filter the diameter data to reduce noise
window = 15;
b = (1/window)*ones(1,window);
a = 1;
dia_avg = zeros(dia_row,dia_col);
for i = 1:dia_row
    dia_avg(i,:) = filtfilt(b,a,dia_mat(i,:));
end

% Initialize the parameters for the wavelet transform (CWT)
dt = 1/15;   % Sampling Frequency (FPS of the video)
dia_t = dt:dt:length(dia_avg)*dt;   % Sampling time points

NumVoices = 64;
a0 = 2^(1/NumVoices);
wavCenterFreq = 5/(2*pi);
minfreq = 0.05;
maxfreq = 1;
minscale = wavCenterFreq/(maxfreq*dt);
maxscale = wavCenterFreq/(minfreq*dt);
minscale = floor(NumVoices*log2(minscale));
maxscale = ceil(NumVoices*log2(maxscale));
scales = a0.^(minscale:maxscale).*dt;

gap = 10;   % Intervals for calculating the 

% Initialize the variables for storing the CWT data
cwt_dia = cell(floor(dia_col/gap),1);
dia_spec = cell(floor(dia_col/gap),1);
abs_dia_spec = cell(floor(dia_col/gap),1);

% Calculate the CWT for the analysis windows
for p = 1:gap:floor(dia_col/gap)*gap
    a = (p-1)/gap+1;
    cwt_dia{a,1} = cwtft({dia_avg(:,p),dt},...
        'wavelet','bump','scales',scales,'padmode','symw');
    dia_spec{a,1} = cwt_dia{a,1}.cfs;   % Calculated Spectrogram
    dia_freq = cwt_dia{a,1}.frequencies;   % Corresponding frequency array
    % Calculate the magnitude of the CWT coefficients
    abs_dia_spec{a,1} = abs(dia_spec{a,1});
end

% Calculate the frequency and phase information from the spectrogram
[row, col] = size(dia_spec{1,1});
freq_max = zeros(col,floor(dia_col/gap));
phase_max = zeros(col,floor(dia_col/gap));
for p = 1:gap:floor(dia_col/gap)*gap
    a = (p-1)/gap+1;
    for i = 1:col
        [dia_freq_peaks, freq] = findpeaks(abs_dia_spec{a,1}(:,i));   % Find the frequencies at which the coefficients are maximum
        [~, freq_loc_max] = max(dia_freq_peaks);
        freq_max(i,a) = dia_freq(freq(freq_loc_max));   % calculate the maximum frequency
        phase_max(i,a) = angle(dia_spec{a,1}(freq(freq_loc_max),i));   % Calculate the phase at the maximum frequency
    end
end

% Save the CWT data
save(strrep(dia_name,'_mat','_avg_cwt'),'dia_avg','dia_spec','abs_dia_spec',...
    'dia_freq','dia_t','freq_max','phase_max','-v7.3')

%% Wavelet Transform for Pressure Difference
% 
% %Read Pressure data
% pressure_file=uigetfile('*.txt');
% pressure_fileID = fopen(pressure_file);
% P = textscan(pressure_fileID,'%f %f %f %f %f %f %f %f %f %f %f %f','EndOfLine','\r\n');
% P = cell2mat(P);
% 
% td=0:1/15:(length(dia_mat)-1)/15;
% td = td';
% tp=P(:,1)-P(1,1);
% actP=zeros(length(td),12);
% 
% actP(1,:)=P(1,:);
% for i=2:length(td)
%     [~, ind]=min(abs(tp-td(i,1)));
%     actP(i,:)=P(ind,:);
% end
% 
% dt = 1/15;   % Sampling Frequency (FPS)
% dia_t = dt:dt:length(dia_avg)*dt;
% 
% NumVoices = 64;
% a0 = 2^(1/NumVoices);
% wavCenterFreq = 5/(2*pi);
% minfreq = 0.05;
% maxfreq = 1;
% minscale = wavCenterFreq/(maxfreq*dt);
% maxscale = wavCenterFreq/(minfreq*dt);
% minscale = floor(NumVoices*log2(minscale));
% maxscale = ceil(NumVoices*log2(maxscale));
% scales = a0.^(minscale:maxscale).*dt;
% 
% cwt_pr = cwtft({actP(:,4),dt},'wavelet','bump','scales',scales,'padmode','symw');
% pr_spec = cwt_pr.cfs;
% pr_freq = cwt_pr.frequencies;
% 
% abs_pr_spec = abs(pr_spec);
% 
% %% Plot Pressure Spectrogram
% 
% f1 = 0.125;
% f2 = 0.06;
% figure
% helperCWTTimeFreqPlot(pr_spec,dia_t,pr_freq,'surf','CWT of Diameter for APSS','Seconds','Hz')
% hold on
% line([0 600],[f1 f1],[10^4 10^4],'Color',[1 1 1],'LineStyle','--')
% line([0 600],[f2 f2],[10^4 10^4],'Color',[1 1 1],'LineStyle','--')
% 
% %% Plot Diameter Spectrogram and Freqeuency vs Time
% 
% f1 = 0.125;
% f2 = 0.06;
% loc = 10;
% figure
% helperCWTTimeFreqPlot(dia_spec{loc,1},dia_t,dia_freq,'surf','CWT of Diameter for APSS','Seconds','Hz')
% hold on
% line([0 600],[f1 f1],[10^4 10^4],'Color',[1 1 1],'LineStyle','--')
% line([0 600],[f2 f2],[10^4 10^4],'Color',[1 1 1],'LineStyle','--')
% 
% figure
% plot(dia_t,freq_max(:,loc),'k-');
% hold on
% line([0 600],[f1 f1],[10^4 10^4],'Color',[0.5 0.5 0.5],'LineStyle','--')
% line([0 600],[f2 f2],[10^4 10^4],'Color',[0.5 0.5 0.5],'LineStyle','--')