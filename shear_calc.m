% Read the saved diameter file
dia_name = uigetfile('*.mat');
load(dia_name);
[dia_row, dia_col] = size(dia_mat);

% Low pass filter the diameter data to reduce noise
window = 20;
b = (1/window)*ones(1,window);
a = 1;
dia_avg = zeros(dia_row,dia_col);
for i = 1:dia_row
    dia_avg(i,:) = filtfilt(b,a,dia_mat(i,:));
end

dt = 1/15;   % Sampling Frequency (FPS of the video)
dia_t = dt:dt:length(dia_avg)*dt;   % Sampling time points
dia_t = dia_t';

%% Calculation of the diastolic diameter

[dia_poly, dia_peak] = dia_envelope(dia_t, dia_avg);   % Peak diameter in pixels detected by polynomial fitting
dpeak = dia_poly*3.4*10^-4;   % Peak diamater in cm
mu = 6.78e-3; % water viscosity @ 38C, dynes-s/cm^2

%% Read the pressure and syringe position data
pressure_file=uigetfile('*.txt');
pressure_fileID = fopen(pressure_file);
P = textscan(pressure_fileID,'%f %f %f %f %f %f %f %f %f %f %f %f','EndOfLine','\r\n');
P = cell2mat(P);

% Adjust the pressure data array to the same size as the diameter array
td=0:1/15:(length(dia_mat)-1)/15;
td = td';
tp=P(:,1)-P(1,1);
actP=zeros(length(td),12);
actP(1,:)=P(1,:);
for i=2:length(td)
    [~, ind]=min(abs(tp-td(i,1)));
    actP(i,:)=P(ind,:);
end

t = actP(:,1);   % Time array
x1 = actP(:,8);   % Position data for syringe 1
x2 = actP(:,9);   % Position data for syringe 2
switched = actP(:,12);   % Position of the solenoid valve

%% Calculation of the wall shear stress

fps = 15;   % FPS of the video
wind_len = 1;   % Length of the calculation window (in sec)
v = ones(length(dia_t)-wind_len*fps,1);   % Initialize the velocity array
Q = ones(length(dia_t)-wind_len*fps,1);   % Initialize the flow rate array
t_diff = ones(length(dia_t)-wind_len*fps,1);   % Initialize the time array
tau = ones(length(dia_t)-wind_len*fps,1);   % Initialize the wall shear stress array
for i = 1:length(dia_t)-wind_len*fps
    x1diff = 0;
    x2diff = 0;
    for j = 1:wind_len*fps
    % Sum up the syringe positions based on the state of the solenoid valve
        if switched(j+i,1) == 0
            x1diff = x1diff + (x1(j+i,1) - x1(j+i-1,1));
            x2diff = x2diff + (x2(j+i-1,1) - x2(j+i,1));
        else
            x1diff = x1diff + (x1(j+i-1,1) - x1(j+i,1));
            x2diff = x2diff + (x2(j+i,1) - x2(j+i-1,1));
        end
    end
    % Average velocity over the specified time window (cm/s)
    v(i,1) = 0.5*(x1diff+x2diff)*10^-1./wind_len;
    % Convert velocity to flow (mL/sec)
    Q(i,1) = 0.25*pi*.146^2*v(i,1);     % inner diameter of tubing is 1.46 mm
    % Place time points in the middle of the calculation windows
    t_diff(i,1) = dia_t((i-1)+floor(wind_len*fps/2),1);
    % Average shear stress over the time window
    tau(i,1) = 4*mu*Q(i,1)./(pi*(dpeak/2).^3);   % Convert mu to dyne-s/cm2
end