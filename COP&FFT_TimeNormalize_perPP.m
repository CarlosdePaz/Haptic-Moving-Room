clc; clear all; 
cd('C:\Investigación\Moving Room\data')
load('ModalityCell.mat')

%% ------------------------- 1.Freq 0.1 COP-Y --------------------------

clc; 
Fs_filter = 50; % sampling frequency
Fc = 7.14; % cut-off frequency
order = 4;
j = 1; n = 1;

for row = 1:45%48

    % ----------------------- COP Vision data -----------------------------
    M_vision = table2array(Vision_01_cell{row}(:,2:2:4));
    M_1_vision = [M_vision(:,1),M_vision(:,end)./1000];
    COP_vision = normalize(M_1_vision(M_1_vision(:,2) > 15 & M_1_vision(:,2) < 110));%,
    [b,a] = butter(order,Fc/(Fs_filter/2),'low');
    COP_cf_vision = filtfilt(b,a,COP_vision);
    COP_cf_n_vision = lengthNormalize(COP_cf_vision, 4750);
    COP__vision(j,:) = COP_cf_n_vision;
    COP_norm_vision = timeNormalize(COP__vision(j,:));
    
    % FFT Vision
    TT_vision = table2array(Vision_01_cell{row}(:,4))./1000;
    T_vision = TT_vision(TT_vision > 5 & TT_vision < 110);        
    tCOP_vision = linspace(T_vision(1),T_vision(end),1000);
    Ts_vision = mean(diff(tCOP_vision));
    Fs_vision = 1/Ts_vision;
    L_vision = length(COP_norm_vision);
    
    % FFT command
    Y = fft(COP_norm_vision);
    P2 = (abs(Y).^2)/(L_vision*Fs_vision);
    P1 = P2(1:L_vision/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs_vision*(0:(L_vision/2))/L_vision;
    
%     f(f > 1) = [];
%     FFT_COP_vision = [f;P1(1:length(f))];
    FFT_COP_vision = [f;P1];
    Freq_FFT_vision = find(FFT_COP_vision(1,:) > .0975 & FFT_COP_vision(1,:) < .105);
    FFT_vision(j,:) = [FFT_COP_vision(1,Freq_FFT_vision), FFT_COP_vision(2,Freq_FFT_vision), median(FFT_COP_vision(2,:))];
    
    clear Y P1 P2 f
    
    % ----------------------- COP Device data -----------------------------
    M_device = table2array(Device_01_cell{row}(:,2:2:4));
    M_1_device = [M_device(:,1),M_device(:,end)./1000];
    COP_device = normalize(M_1_device(M_1_device(:,2) > 15 & M_1_device(:,2) < 110));%,
    [b,a] = butter(order,Fc/(Fs_filter/2),'low');
    COP_cf_device = filtfilt(b,a,COP_device);
    COP_cf_n_device = lengthNormalize(COP_cf_device, 4750);
    COP__device(j,:) = COP_cf_n_device;
    COP_norm_device = timeNormalize(COP__device(j,:));
    
    % FFT Device
    TT_device = table2array(Device_01_cell{row}(:,4))./1000;
    T_device = TT_device(TT_device > 5 & TT_device < 110);        
    tCOP_device = linspace(T_device(1),T_device(end),1000);
    Ts_device = mean(diff(tCOP_device));
    Fs_device = 1/Ts_device;
    L_device = length(COP_norm_device);
    
    % FFT command
    Y = fft(COP_norm_device);
    P2 = (abs(Y).^2)/(L_device*Fs_device);
    P1 = P2(1:L_device/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs_device*(0:(L_device/2))/L_device;
    
%     f(f > 1) = [];
%     FFT_COP_device = [f;P1(1:length(f))];

    FFT_COP_device = [f;P1];
    Freq_FFT_device = find(FFT_COP_device(1,:) > .0975 & FFT_COP_device(1,:) < .105);
    FFT_device(j,:) = [FFT_COP_device(1,Freq_FFT_device), FFT_COP_device(2,Freq_FFT_device), median(FFT_COP_device(2,:))];
    
    clear Y P1 P2 f

    % ----------------------- COP Combination data -----------------------------
    M_combination = table2array(Combination_01_cell{row}(:,2:2:4));
    M_1_combination = [M_combination(:,1),M_combination(:,end)./1000];
    COP_combination = normalize(M_1_combination(M_1_combination(:,2) > 15 & M_1_combination(:,2) < 110));%,
    [b,a] = butter(order,Fc/(Fs_filter/2),'low');
    COP_cf_combination = filtfilt(b,a,COP_combination);
    COP_cf_n_combination = lengthNormalize(COP_cf_combination, 4750);
    COP__combination(j,:) = COP_cf_n_combination;
    COP_norm_combination = timeNormalize(COP__combination(j,:));
    
    % FFT combination
    TT_combination = table2array(Combination_01_cell{row}(:,4))./1000;
    T_combination = TT_combination(TT_combination > 5 & TT_combination < 110);        
    tCOP_combination = linspace(T_combination(1),T_combination(end),1000);
    Ts_combination = mean(diff(tCOP_combination));
    Fs_combination = 1/Ts_combination;
    L_combination = length(COP_norm_combination);
    
    % FFT command
    Y = fft(COP_norm_combination);
    P2 = (abs(Y).^2)/(L_combination*Fs_combination);
    P1 = P2(1:L_combination/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs_combination*(0:(L_combination/2))/L_combination;
    
%     f(f > 1) = [];
%     FFT_COP_combination = [f;P1(1:length(f))];

    FFT_COP_combination = [f;P1];
    Freq_FFT_combination = find(FFT_COP_combination(1,:) > .0975 & FFT_COP_combination(1,:) < .105);
    FFT_combination(j,:) = [FFT_COP_combination(1,Freq_FFT_combination), FFT_COP_combination(2,Freq_FFT_combination), median(FFT_COP_combination(2,:))];

    j = j+1;
    if j > 3
        COP_avg_vision(n,:) = mean(COP__vision(1:3,:));
        FFT_avg_vision(n,:) = mean(FFT_vision(1:3,:));
        COP_avg_device(n,:) = mean(COP__device(1:3,:));
        FFT_avg_device(n,:) = mean(FFT_device(1:3,:));
        COP_avg_combination(n,:) = mean(COP__combination(1:3,:));
        FFT_avg_combination(n,:) = mean(FFT_combination(1:3,:));
        j = 1;
        n = n+1;
    end
    
    clearvars -except j Fs_filter Fc order n FFT_vision Vision_01_cell FFT_avg_vision COP_avg_vision FFT_device Device_01_cell FFT_avg_device COP_avg_device FFT_combination Combination_01_cell FFT_avg_combination COP_avg_combination
end

COP_avg_vision(6,:) = [];
COP_avg_device(6,:) = [];
COP_avg_combination(6,:) = [];

FFT_avg_vision(6,:) = [];
FFT_avg_device(6,:) = [];
FFT_avg_combination(6,:) = [];

FFT_01 = [FFT_avg_vision, FFT_avg_device, FFT_avg_combination];
FFT_Table_01 = array2table(FFT_01, 'VariableNames', {'F01_V','Power01_V','MedianPower01_V','F01_D','Power01_D','MedianPower01_D','F01_C','Power01_C','MedianPower01_C'});

%%
clc; clear all; 
cd('C:\Investigación\Moving Room\data')
load('ModalityCells.mat')

%% ------------------------- 2.Freq 0.75 COP-Y --------------------------

clc; 
Fs_filter = 50; % sampling frequency
Fc = 7.14; % cut-off frequency
order = 4;
j = 1; n = 1;
[b,a] = butter(order,Fc/(Fs_filter/2),'low');

for row = 1:48

    % ----------------------- COP Vision data -----------------------------
    M_vision = table2array(Vision_075_cell{row}(:,2:2:4));
    M_1_vision = [M_vision(:,1),M_vision(:,end)./1000];
    COP_vision = normalize(M_1_vision(M_1_vision(:,2) > (5+1.34) & M_1_vision(:,2) < (M_1_vision(end,2)-5-1.34*1.5)));%,'center'
    COP_cf_vision = filtfilt(b,a,COP_vision);
    COP_cf_n_vision = lengthNormalize(COP_cf_vision, 640);
    COP__vision(j,:) = COP_cf_n_vision;
    COP_norm_vision = timeNormalize(COP_avg_vision(j,:));
    
    % FFT Vision
    TT_vision = table2array(Vision_075_cell{row}(:,4))./1000;
    T_vision = TT_vision(TT_vision > (5+1.34) & TT_vision < (TT_vision(end)-5-1.34*1.5));        
    tCOP_vision = linspace(T_vision(1),T_vision(end),1000);
    Ts_vision = mean(diff(tCOP_vision));
    Fs_vision = 1/Ts_vision;
    L_vision = length(COP_norm_vision);
    
    % FFT command
    Y = fft(COP_norm_vision);
    P2 = (abs(Y).^2)/(L_vision*Fs_vision);
    P1 = P2(1:L_vision/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs_vision*(0:(L_vision/2))/L_vision;
    
%     f(f > 1) = [];
%     FFT_COP_vision = [f;P1(1:length(f))];
    FFT_COP_vision = [f;P1];
    Freq_FFT_vision = find(FFT_COP_vision(1,:) > .70 & FFT_COP_vision(1,:) < .76);
    FFT_vision(j,:) = [FFT_COP_vision(1,Freq_FFT_vision), FFT_COP_vision(2,Freq_FFT_vision), median(FFT_COP_vision(2,:))];
    
    clear Y P1 P2 f
    
    % ----------------------- COP Device data -----------------------------
    M_device = table2array(Device_075_cell{row}(:,2:2:4));
    M_1_device = [M_device(:,1),M_device(:,end)./1000];
    COP_device = normalize(M_1_device(M_1_device(:,2) > (5+1.34) & M_1_device(:,2) < (M_1_device(end,2)-5-1.34*1.5)));%,'center'
    COP_cf_device = filtfilt(b,a,COP_device);
    COP_cf_n_device = lengthNormalize(COP_cf_device, 640);
    COP__device(j,:) = COP_cf_n_device;
    COP_norm_device = timeNormalize(COP_avg_device(j,:));
    
    % FFT Device
    TT_device = table2array(Device_075_cell{row}(:,4))./1000;
    T_device = TT_device(TT_device > (5+1.34) & TT_device < (TT_device(end)-5-1.34*1.5));        
    tCOP_device = linspace(T_device(1),T_device(end),1000);
    Ts_device = mean(diff(tCOP_device));
    Fs_device = 1/Ts_device;
    L_device = length(COP_norm_device);
    
    % FFT command
    Y = fft(COP_norm_device);
    P2 = (abs(Y).^2)/(L_device*Fs_device);
    P1 = P2(1:L_device/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs_device*(0:(L_device/2))/L_device;
    
%     f(f > 1) = [];
%     FFT_COP_device = [f;P1(1:length(f))];
    FFT_COP_device = [f;P1];
    Freq_FFT_device = find(FFT_COP_device(1,:) > .70 & FFT_COP_device(1,:) < .76);
    FFT_device(j,:) = [FFT_COP_device(1,Freq_FFT_device), FFT_COP_device(2,Freq_FFT_device), median(FFT_COP_device(2,:))];

    clear Y P1 P2 f
    % ----------------------- COP Combination data -----------------------------
    M_combination = table2array(Combination_075_cell{row}(:,2:2:4));
    M_1_combination = [M_combination(:,1),M_combination(:,end)./1000];
    COP_combination = normalize(M_1_combination(M_1_combination(:,2) > (5+1.34) & M_1_combination(:,2) < (M_1_combination(end,2)-5-1.34*1.5)));%,'center'
    COP_cf_combination = filtfilt(b,a,COP_combination);
    COP_cf_n_combination = lengthNormalize(COP_cf_combination, 640);
    COP__combination(j,:) = COP_cf_n_combination;
    COP_norm_combination = timeNormalize(COP_avg_combination(j,:));
    
    % FFT combination
    TT_combination = table2array(Combination_075_cell{row}(:,4))./1000;
    T_combination = TT_combination(TT_combination > (5+1.34) & TT_combination < (TT_combination(end)-5-1.34*1.5));        
    tCOP_combination = linspace(T_combination(1),T_combination(end),1000);
    Ts_combination = mean(diff(tCOP_combination));
    Fs_combination = 1/Ts_combination;
    L_combination = length(COP_norm_combination);
    
    % FFT command
    Y = fft(COP_norm_combination);
    P2 = (abs(Y).^2)/(L_combination*Fs_combination);
    P1 = P2(1:L_combination/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs_combination*(0:(L_combination/2))/L_combination;
    
%     f(f > 1) = [];
%     FFT_COP_combination = [f;P1(1:length(f))];
    FFT_COP_combination = [f;P1];
    Freq_FFT_combination = find(FFT_COP_combination(1,:) > .70 & FFT_COP_combination(1,:) < .76);
    FFT_combination(j,:) = [FFT_COP_combination(1,Freq_FFT_combination), FFT_COP_combination(2,Freq_FFT_combination), median(FFT_COP_combination(2,:))];

    j = j+1;
    if j > 3
        COP_avg_vision(n,:) = mean(COP__vision(1:3,:));
        FFT_avg_vision(n,:) = mean(FFT_vision(1:3,:));
        COP_avg_device(n,:) = mean(COP__device(1:3,:));
        FFT_avg_device(n,:) = mean(FFT_device(1:3,:));
        COP_avg_combination(n,:) = mean(COP__combination(1:3,:));
        FFT_avg_combination(n,:) = mean(FFT_combination(1:3,:));
        j = 1;
        n = n+1;
    end
    
    clearvars -except j a b Fs_filter Fc order n FFT_vision Vision_075_cell FFT_avg_vision COP_avg_vision FFT_device Device_075_cell FFT_avg_device COP_avg_device FFT_combination Combination_075_cell FFT_avg_combination COP_avg_combination
end

COP_avg_vision(6,:) = [];
COP_avg_device(6,:) = [];
COP_avg_combination(6,:) = [];

FFT_avg_vision(6,:) = [];
FFT_avg_device(6,:) = [];
FFT_avg_combination(6,:) = [];

FFT_075 = [FFT_avg_vision, FFT_avg_device, FFT_avg_combination];
FFT_Table_075 = array2table(FFT_075, 'VariableNames', {'F075_V','Power075_V','MedianPower075_V','F075_D','Power075_D','MedianPower075_D','F075_C','Power075_C','MedianPower075_C'});


%%
clc; clear all; 
cd('C:\Investigación\Moving Room\data')
load('ModalityCells.mat')

%% ------------------------- 3. Still COP-Y --------------------------

clc; 
Fs_filter = 50; % sampling frequency
Fc = 7.14; % cut-off frequency
order = 4;
j = 1; n = 1;
    [b,a] = butter(order,Fc/(Fs_filter/2),'low');

for row = 1:16

    % ----------------------- COP Vision data -----------------------------
    M_vision = table2array(Vision_still_cell{row}(:,2:2:4));
    M_1_vision = [M_vision(:,1),M_vision(:,end)./1000];
    COP_vision = (M_1_vision(M_1_vision(:,2) > 5 & M_1_vision(:,2) < 35));%,'center'
    COP_cf_vision = filtfilt(b,a,COP_vision);
    COP_cf_n_vision = lengthNormalize(COP_cf_vision, 1500);
    COP_avg_vision(row,:) = COP_cf_n_vision;
    COP_norm_vision = timeNormalize(COP_avg_vision(row,:));
    
    % FFT Vision
    TT_vision = table2array(Vision_still_cell{row}(:,4))./1000;
    T_vision = TT_vision(TT_vision > 5 & TT_vision < 35);        
    tCOP_vision = linspace(T_vision(1),T_vision(end),1000);
    Ts_vision = mean(diff(tCOP_vision));
    Fs_vision = 1/Ts_vision;
    L_vision = length(COP_norm_vision);
    
    % FFT command
    Y = fft(COP_norm_vision);
    P2 = (abs(Y).^2)/(L_vision*Fs_vision);
    P1 = P2(1:L_vision/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs_vision*(0:(L_vision/2))/L_vision;
    
%     f(f > 1) = [];
%     FFT_COP_vision = [f;P1(1:length(f))];
    FFT_COP_vision = [f;P1];
    Freq_FFT_vision_01 = find(FFT_COP_vision(1,:) > .0975 & FFT_COP_vision(1,:) < .105);
    Freq_FFT_vision_075 = find(FFT_COP_vision(1,:) > .70 & FFT_COP_vision(1,:) < .76);

    FFT_avg_vision(row,:) = [FFT_COP_vision(1,Freq_FFT_vision_01), FFT_COP_vision(2,Freq_FFT_vision_01), FFT_COP_vision(1,Freq_FFT_vision_075), FFT_COP_vision(2,Freq_FFT_vision_075),median(FFT_COP_vision(2,:))];
    
    clear Y P1 P2 f
    % ----------------------- COP Device data -----------------------------
    M_device = table2array(Device_still_cell{row}(:,2:2:4));
    M_1_device = [M_device(:,1),M_device(:,end)./1000];
    COP_device = (M_1_device(M_1_device(:,2) > 5 & M_1_device(:,2) < 35));%,'center'
    COP_cf_device = filtfilt(b,a,COP_device);
    COP_cf_n_device = lengthNormalize(COP_cf_device, 1500);
    COP_avg_device(row,:) = COP_cf_n_device;
    COP_norm_device = timeNormalize(COP_avg_device(row,:));
    
    % FFT Device
    TT_device = table2array(Device_still_cell{row}(:,4))./1000;
    T_device = TT_device(TT_device > 5 & TT_device < 35);        
    tCOP_device = linspace(T_device(1),T_device(end),1000);
    Ts_device = mean(diff(tCOP_device));
    Fs_device = 1/Ts_device;
    L_device = length(COP_norm_device);
    
    % FFT command
    Y = fft(COP_norm_device);
    P2 = (abs(Y).^2)/(L_device*Fs_device);
    P1 = P2(1:L_device/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs_device*(0:(L_device/2))/L_device;
    
%     f(f > 1) = [];
%     FFT_COP_device = [f;P1(1:length(f))];
    FFT_COP_device = [f;P1];
    Freq_FFT_device_01 = find(FFT_COP_device(1,:) > .0975 & FFT_COP_device(1,:) < .105);
    Freq_FFT_device_075 = find(FFT_COP_device(1,:) > .70 & FFT_COP_device(1,:) < .76);

    FFT_avg_device(row,:) = [FFT_COP_device(1,Freq_FFT_device_01), FFT_COP_device(2,Freq_FFT_device_01), FFT_COP_device(1,Freq_FFT_device_075), FFT_COP_device(2,Freq_FFT_device_075),median(FFT_COP_device(2,:))];
    
    clear Y P1 P2 f
    % ----------------------- COP Combination data -----------------------------
    M_combination = table2array(Combination_still_cell{row}(:,2:2:4));
    M_1_combination = [M_combination(:,1),M_combination(:,end)./1000];
    COP_combination = (M_1_combination(M_1_combination(:,2) > 5 & M_1_combination(:,2) < 35));%,'center'
    [b,a] = butter(order,Fc/(Fs_filter/2),'low');
    COP_cf_combination = filtfilt(b,a,COP_combination);
    COP_cf_n_combination = lengthNormalize(COP_cf_combination, 1500);
    COP_avg_combination(row,:) = COP_cf_n_combination;
    COP_norm_combination = timeNormalize(COP_avg_combination(row,:));
    
    % FFT combination
    TT_combination = table2array(Combination_still_cell{row}(:,4))./1000;
    T_combination = TT_combination(TT_combination > 5 & TT_combination < 35);        
    tCOP_combination = linspace(T_combination(1),T_combination(end),1000);
    Ts_combination = mean(diff(tCOP_combination));
    Fs_combination = 1/Ts_combination;
    L_combination = length(COP_norm_combination);
    
    % FFT command
    Y = fft(COP_norm_combination);
    P2 = (abs(Y).^2)/(L_combination*Fs_combination);
    P1 = P2(1:L_combination/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs_combination*(0:(L_combination/2))/L_combination;
    
%     f(f > 1) = [];
    FFT_COP_combination = [f;P1];
    Freq_FFT_combination_01 = find(FFT_COP_combination(1,:) > .0975 & FFT_COP_combination(1,:) < .105);
    Freq_FFT_combination_075 = find(FFT_COP_combination(1,:) > .70 & FFT_COP_combination(1,:) < .76);
    FFT_avg_combination(row,:) = [FFT_COP_combination(1,Freq_FFT_combination_01), FFT_COP_combination(2,Freq_FFT_combination_01), FFT_COP_combination(1,Freq_FFT_combination_075), FFT_COP_combination(2,Freq_FFT_combination_075),median(FFT_COP_combination(2,:))];

    clear Y P1 P2 f
    
    clearvars -except b a j Fs_filter Fc order n FFT_vision Vision_still_cell FFT_avg_vision COP_avg_vision FFT_device Device_still_cell FFT_avg_device COP_avg_device FFT_combination Combination_still_cell FFT_avg_combination COP_avg_combination
end

COP_avg_vision(6,:) = [];
COP_avg_device(6,:) = [];
COP_avg_combination(6,:) = [];

FFT_avg_vision(6,:) = [];
FFT_avg_device(6,:) = [];
FFT_avg_combination(6,:) = [];

FFT_still = [FFT_avg_vision, FFT_avg_device, FFT_avg_combination];
FFT_Table_still = array2table(FFT_still, 'VariableNames', {'F01_Still_V','Power01_Still_V','F075_Still_V','Power075_Still_V','MedianPowerStill_V', 'F01_Still_D','Power01_Still_D','F075_Still_D','Power075_Still_D','MedianPowerStill_D', 'F01_Still_C','Power01_Still_C','F075_Still_C','Power075_Still_C','MedianPowerStill_C'});

%%

FFT_COP = [FFT_01, FFT_075, FFT_still];

names = {'F01_V','Power01_V','MedianPower01_V','F01_D','Power01_D','MedianPower01_D','F01_C','Power01_C','MedianPower01_C','F075_V','Power075_V','MedianPower075_V','F075_D','Power075_D','MedianPower075_D','F075_C','Power075_C','MedianPower075_C','F01_Still_V','Power01_Still_V','F075_Still_V','Power075_Still_V','MedianPowerStill_V', 'F01_Still_D','Power01_Still_D','F075_Still_D','Power075_Still_D','MedianPowerStill_D', 'F01_Still_C','Power01_Still_C','F075_Still_C','Power075_Still_C','MedianPowerStill_C'};

FFT_COP_Table = array2table(FFT_COP, 'VariableNames', names);

%% 

FFT_W = [FFT_COP(1:9,2:3),FFT_COP(1:9,5:6),FFT_COP(1:9,8:9),FFT_COP(1:9,11:12),FFT_COP(1:9,14:15),FFT_COP(1:9,17:18)]
FFT_W_mean = [mean(FFT_COP(1:9,2:3)),mean(FFT_COP(1:9,5:6)),mean(FFT_COP(1:9,8:9)),mean(FFT_COP(1:9,11:12)),mean(FFT_COP(1:9,14:15)),mean(FFT_COP(1:9,17:18))];
n = 1;
for i = 1:2:12
    p(n,1) = ranksum(FFT_W(:,i),FFT_W(:,i+1));
    n = n+1;
end


%% Vision

p_Val(1,1) = ranksum(FFT_COP(1:9,2),FFT_COP(1:9,3)); % 1. Power 01 vs Median Power 01
p_Val(2,1) = ranksum(FFT_COP(1:9,11),FFT_COP(1:9,12)); % 2. Power 075 vs Median Power 075 
p_Val(3,1) = ranksum(FFT_COP(1:9,2),FFT_COP(1:9,23)); % 3. Power 01 vs Still Median Power 075 
p_Val(4,1) = ranksum(FFT_COP(1:9,11),FFT_COP(1:9,23)); % 4. Power 075 vs Still Median Power 075 
p_Val(5,1) = ranksum(FFT_COP(1:9,2),FFT_COP(1:9,20)); % 5. Power 01 vs Still 01
p_Val(6,1) = ranksum(FFT_COP(1:9,11),FFT_COP(1:9,22)); % 5. Power 075 vs Still 01

% Device
p_Val(1,2) = ranksum(FFT_COP(1:9,5),FFT_COP(1:9,6)); % 1. Power 01 vs Median Power 01
p_Val(2,2) = ranksum(FFT_COP(1:9,14),FFT_COP(1:9,15)); % 2. Power 075 vs Median Power 075 
p_Val(3,2) = ranksum(FFT_COP(1:9,5),FFT_COP(1:9,28)); % 3. Power 01 vs Still Median Power 075 
p_Val(4,2) = ranksum(FFT_COP(1:9,14),FFT_COP(1:9,28)); % 4. Power 075 vs Still Median Power 075 
p_Val(5,1) = ranksum(FFT_COP(1:9,5),FFT_COP(1:9,25)); % 5. Power 01 vs Still 01
p_Val(6,1) = ranksum(FFT_COP(1:9,14),FFT_COP(1:9,27)); % 5. Power 075 vs Still 01

% Combination
p_Val(1,3) = ranksum(FFT_COP(1:9,5),FFT_COP(1:9,6)); % 1. Power 01 vs Median Power 01
p_Val(2,3) = ranksum(FFT_COP(1:9,14),FFT_COP(1:9,15)); % 2. Power 075 vs Median Power 075 
p_Val(3,3) = ranksum(FFT_COP(1:9,5),FFT_COP(1:9,28)); % 3. Power 01 vs Still Median Power 075 
p_Val(4,3) = ranksum(FFT_COP(1:9,14),FFT_COP(1:9,28)); % 4. Power 075 vs Still Median Power 075 
p_Val(5,3) = ranksum(FFT_COP(1:9,5),FFT_COP(1:9,25)); % 5. Power 01 vs Still 01
p_Val(6,3) = ranksum(FFT_COP(1:9,14),FFT_COP(1:9,27)); % 5. Power 075 vs Still 01

%% Vision
p_Val(1,1) = ranksum(FFT_COP_9(:,2),FFT_COP_9(:,3)); % 1. Power 01 vs Median Power 01
p_Val(2,1) = ranksum(FFT_COP_9(:,11),FFT_COP_9(:,12)); % 2. Power 075 vs Median Power 075 
p_Val(3,1) = ranksum(FFT_COP_9(:,2),FFT_COP_9(:,23)); % 3. Power 01 vs Still Median Power 075 
p_Val(4,1) = ranksum(FFT_COP_9(:,11),FFT_COP_9(:,23)); % 4. Power 075 vs Still Median Power 075 
p_Val(5,1) = ranksum(FFT_COP_9(:,2),FFT_COP_9(:,20)); % 5. Power 01 vs Still 01
p_Val(6,1) = ranksum(FFT_COP_9(:,11),FFT_COP_9(:,22)); % 5. Power 075 vs Still 075

% Device
p_Val(1,2) = ranksum(FFT_COP_9(:,5),FFT_COP_9(:,6)); % 1. Power 01 vs Median Power 01
p_Val(2,2) = ranksum(FFT_COP_9(:,14),FFT_COP_9(:,15)); % 2. Power 075 vs Median Power 075 
p_Val(3,2) = ranksum(FFT_COP_9(:,5),FFT_COP_9(:,28)); % 3. Power 01 vs Still Median Power 075 
p_Val(4,2) = ranksum(FFT_COP_9(:,14),FFT_COP_9(:,28)); % 4. Power 075 vs Still Median Power 075 
p_Val(5,1) = ranksum(FFT_COP_9(:,5),FFT_COP_9(:,25)); % 5. Power 01 vs Still 01
p_Val(6,1) = ranksum(FFT_COP_9(:,14),FFT_COP_9(:,27)); % 5. Power 075 vs Still 075

% Combination
p_Val(1,3) = ranksum(FFT_COP_9(:,5),FFT_COP_9(:,6)); % 1. Power 01 vs Median Power 01
p_Val(2,3) = ranksum(FFT_COP_9(:,14),FFT_COP_9(:,15)); % 2. Power 075 vs Median Power 075 
p_Val(3,3) = ranksum(FFT_COP_9(:,5),FFT_COP_9(:,28)); % 3. Power 01 vs Still Median Power 075 
p_Val(4,3) = ranksum(FFT_COP_9(:,14),FFT_COP_9(:,28)); % 4. Power 075 vs Still Median Power 075 
p_Val(5,3) = ranksum(FFT_COP_9(:,5),FFT_COP_9(:,25)); % 5. Power 01 vs Still 01
p_Val(6,3) = ranksum(FFT_COP_9(:,14),FFT_COP_9(:,27)); % 5. Power 075 vs Still 075

