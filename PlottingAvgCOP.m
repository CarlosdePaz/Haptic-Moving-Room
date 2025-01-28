%% Average Plotting

%% Vision 0.1 Hz

clc; clear all; close all;
cd('C:\Investigación\Moving Room\data')
load('ModalityCells.mat')

Fs = 50; % sampling frequency
Fc = 7.14; % cut-off frequency
order = 4;
j = 1; n = 1;

subplot(3,1,1)
figure(1)
hold on
box on
axis([0,1000, -1, 2])
for row = 1:48
    M = table2array(Vision_01_cell{row}(:,2:2:4));
    M_1 = [M(:,1),M(:,end)./1000];
    M_3 = M_1(M_1(:,2) > 15 & M_1(:,2) < 110);
    COP = timeNormalize(M_3);

    % Get coefficients of the filter
    [b,a] = butter(order,Fc/(Fs/2),'low');
    COP_cf = filtfilt(b,a,COP);
    COP_n = normalize(COP_cf,'center');
    plot(COP_n, 'r')
    COP_avg(j,:) = COP_n;
    
    M_Room = [table2array(Vision_Room_01_cell{row}(:,4:4:8))];
    M_1_Room = [M_Room(:,1),M_Room(:,end)./1000];
    M_3_Room = M_1_Room(M_1_Room(:,2) > 15 & M_1_Room(:,2) < 110);
    Room = normalize(timeNormalize(M_3_Room));
    plot(Room, 'b')
    Room_avg(j,:) = Room;
    
    clear M M_1 M_3 COP COP_cf COP_n
    j = j+1;
    if j > 3
        COP_avg(n,:) = mean(COP_avg(1:3,:));
        Room_avg(n,:) = mean(Room_avg(1:3,:));
        j = 1;
        n = n+1;
    end
end
hold off

COP_avg(6,:) = [];

subplot(3,1,2)
hold on
box on
for pp = 1:15
    plot(COP_avg(pp,:), 'r')
    plot(Room_avg(pp,:), 'b')
end
hold off

subplot(3,1,3)
hold on
box on
COP_avg_total = mean(COP_avg);
Room_avg_total = mean(Room_avg);
plot(COP_avg_total,'r')
plot(Room_avg_total, 'b')
hold off

COP_V01 = COP_avg;
COP_V01_total = COP_avg_total;
Room_V01 = Room_avg;
Room_V01_total = Room_avg_total;

clearvars -except COP_V01 COP_V01_total Room_V01 Room_V01_total Vision_01_cell Device_01_cell Combination_01_cell Vision_Room_01_cell Device_Room_01_cell Combination_Room_01_cell

% Device 0.1 Hz

cd('C:\Investigación\Moving Room\data')

Fs = 50; % sampling frequency
Fc = 7.14; % cut-off frequency
order = 4;
j = 1; n = 1;

subplot(3,1,1)
figure(1)
hold on
box on
axis([0,1000, -1, 2])
for row = 1:48
    M = table2array(Device_01_cell{row}(:,2:2:4));
    M_1 = [M(:,1),M(:,end)./1000];
    M_3 = M_1(M_1(:,2) > 15 & M_1(:,2) < 110);
    COP = timeNormalize(M_3);

    % Get coefficients of the filter
    [b,a] = butter(order,Fc/(Fs/2),'low');
    COP_cf = filtfilt(b,a,COP);
    COP_n = normalize(COP_cf,'center');
    plot(COP_n, 'r')
    COP_avg(j,:) = COP_n;
    
    M_Room = [table2array(Device_Room_01_cell{row}(:,4:4:8))];
    M_1_Room = [M_Room(:,1),M_Room(:,end)./1000];
    M_3_Room = M_1_Room(M_1_Room(:,2) > 15 & M_1_Room(:,2) < 110);
    Room = normalize(timeNormalize(M_3_Room));
    plot(Room, 'b')
    Room_avg(j,:) = Room;
    
    clear M M_1 M_3 COP COP_cf COP_n
    j = j+1;
    if j > 3
        COP_avg(n,:) = mean(COP_avg(1:3,:));
        Room_avg(n,:) = mean(Room_avg(1:3,:));
        j = 1;
        n = n+1;
    end
end
hold off

COP_avg(6,:) = [];

subplot(3,1,2)
hold on
box on
for pp = 1:15
    plot(COP_avg(pp,:), 'r')
    plot(Room_avg(pp,:), 'b')
end
hold off

subplot(3,1,3)
hold on
box on
COP_avg_total = mean(COP_avg);
Room_avg_total = mean(Room_avg);
plot(COP_avg_total,'r')
plot(Room_avg_total, 'b')
hold off

COP_D01 = COP_avg;
COP_D01_total = COP_avg_total;
Room_D01 = Room_avg;
Room_D01_total = Room_avg_total;

clearvars -except COP_V01 COP_V01_total Room_V01 Room_V01_total COP_D01 COP_D01_total Room_D01 Room_D01_total Vision_01_cell Device_01_cell Combination_01_cell Vision_Room_01_cell Device_Room_01_cell Combination_Room_01_cell

% Combination 0.1 Hz

cd('C:\Investigación\Moving Room\data')

Fs = 50; % sampling frequency
Fc = 7.14; % cut-off frequency
order = 4;
j = 1; n = 1;

subplot(3,1,1)
figure(1)
hold on
box on
axis([0,1000, -1, 2])
for row = 1:48
    M = table2array(Combination_01_cell{row}(:,2:2:4));
    M_1 = [M(:,1),M(:,end)./1000];
    M_3 = M_1(M_1(:,2) > 15 & M_1(:,2) < 110);
    COP = timeNormalize(M_3);

    % Get coefficients of the filter
    [b,a] = butter(order,Fc/(Fs/2),'low');
    COP_cf = filtfilt(b,a,COP);
    COP_n = normalize(COP_cf,'center');
    plot(COP_n, 'r')
    COP_avg(j,:) = COP_n;
    
    M_Room = [table2array(Combination_Room_01_cell{row}(:,4:4:8))];
    M_1_Room = [M_Room(:,1),M_Room(:,end)./1000];
    M_3_Room = M_1_Room(M_1_Room(:,2) > 15 & M_1_Room(:,2) < 110);
    Room = normalize(timeNormalize(M_3_Room));
    plot(Room, 'b')
    Room_avg(j,:) = Room;
    
    clear M M_1 M_3 COP COP_cf COP_n
    j = j+1;
    if j > 3
        COP_avg(n,:) = mean(COP_avg(1:3,:));
        Room_avg(n,:) = mean(Room_avg(1:3,:));
        j = 1;
        n = n+1;
    end
end
hold off

COP_avg(6,:) = [];

subplot(3,1,2)
hold on
box on
for pp = 1:15
    plot(COP_avg(pp,:), 'r')
    plot(Room_avg(pp,:), 'b')
end
hold off

subplot(3,1,3)
hold on
box on
COP_avg_total = mean(COP_avg);
Room_avg_total = mean(Room_avg);
plot(COP_avg_total,'r')
plot(Room_avg_total, 'b')
hold off

COP_C01 = COP_avg;
COP_C01_total = COP_avg_total;
Room_C01 = Room_avg;
Room_C01_total = Room_avg_total;

Room_01_total = mean([Room_V01_total;Room_D01_total;Room_C01_total])./10;

clearvars -except COP_V01 COP_V01_total Room_V01 Room_V01_total COP_D01 COP_D01_total Room_D01 Room_D01_total COP_C01 COP_C01_total Room_C01 Room_C01_total Room_01_total Vision_01_cell Device_01_cell Combination_01_cell Vision_Room_01_cell Device_Room_01_cell Combination_Room_01_cell
%%
close all
figure()
hold on
box on
title('COP 0.1 Hz')
plot(COP_V01_total, 'b', 'LineWidth', 1.5)
plot(COP_D01_total, 'r','LineWidth', 1.5)
plot(COP_C01_total, 'g','LineWidth', 1.5)
plot(Room_01_total, 'k','LineWidth', 1.5)
legend('Vision', 'Device', 'Combination','Wall', 'Location', 'SouthEast')
hold off



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Vision 0.75 Hz

clc; clear all; close all;
cd('C:\Investigación\Moving Room\data')
load('ModalityCells.mat')

Fs = 50; % sampling frequency
Fc = 7.14; % cut-off frequency
order = 4;
j = 1; n = 1;

subplot(3,1,1)
figure(1)
hold on
box on
axis([0,1000, -1, 2])
for row = 1:48
    M = table2array(Vision_075_cell{row}(:,2:2:4));
    M_1 = [M(:,1),M(:,end)./1000];
    M_3 = M_1(M_1(:,2) > (5+1.34) & M_1(:,2) < (M_1(end,2)-5-1.34*1.5));
    COP = timeNormalize(M_3);

    % Get coefficients of the filter
    [b,a] = butter(order,Fc/(Fs/2),'low');
    COP_cf = filtfilt(b,a,COP);
    COP_n = normalize(COP_cf,'center');
    plot(COP_n, 'b')
    COP_avg(j,:) = COP_n;
    
    M_Room = [table2array(Vision_Room_075_cell{row}(:,4:4:8))];
    M_1_Room = [M_Room(:,1),M_Room(:,end)./1000];
    M_3_Room = M_1_Room(M_1_Room(:,2) > (5+1.34) & M_1_Room(:,2) < (M_1_Room(end,2)-5-1.34*1.5));
    Room = normalize(timeNormalize(M_3_Room));
    plot(Room, 'b--')
    Room_avg(j,:) = Room;
    
    clear M M_1 M_3 COP COP_cf COP_n
    j = j+1;
    if j > 3
        COP_avg(n,:) = mean(COP_avg(1:3,:));
        Room_avg(n,:) = mean(Room_avg(1:3,:));
        j = 1;
        n = n+1;
    end
end
hold off

COP_avg(6,:) = [];

subplot(3,1,2)
hold on
box on
for pp = 1:15
    plot(COP_avg(pp,:), 'b')
    plot(Room_avg(pp,:), 'b--')
end
hold off

subplot(3,1,3)
hold on
box on
COP_avg_total = mean(COP_avg);
Room_avg_total = mean(Room_avg);
plot(COP_avg_total,'b')
plot(Room_avg_total, 'b--')
hold off

COP_V075 = COP_avg;
COP_V075_total = COP_avg_total;
Room_V075 = Room_avg;
Room_V075_total = Room_avg_total;

clearvars -except COP_V075 COP_V075_total Room_V075 Room_V075_total

% Device 0.75 Hz

cd('C:\Investigación\Moving Room\data')
load('ModalityCells.mat')

Fs = 50; % sampling frequency
Fc = 7.14; % cut-off frequency
order = 4;
j = 1; n = 1;

subplot(3,1,1)
figure(1)
hold on
box on
axis([0,1000, -1, 2])
for row = 1:48
    M = table2array(Device_075_cell{row}(:,2:2:4));
    M_1 = [M(:,1),M(:,end)./1000];
    M_3 = M_1(M_1(:,2) > (5+1.34) & M_1(:,2) < (M_1(end,2)-5-1.34*1.5));
    COP = timeNormalize(M_3);

    % Get coefficients of the filter
    [b,a] = butter(order,Fc/(Fs/2),'low');
    COP_cf = filtfilt(b,a,COP);
    COP_n = normalize(COP_cf,'center');
    plot(COP_n, 'r')
    COP_avg(j,:) = COP_n;
    
    M_Room = [table2array(Device_Room_075_cell{row}(:,4:4:8))];
    M_1_Room = [M_Room(:,1),M_Room(:,end)./1000];
    M_3_Room = M_1_Room(M_1_Room(:,2) > (5+1.34) & M_1_Room(:,2) < (M_1_Room(end,2)-5-1.34*1.5));
    Room = normalize(timeNormalize(M_3_Room));
    plot(Room, 'r--')
    Room_avg(j,:) = Room;
    
    clear M M_1 M_3 COP COP_cf COP_n
    j = j+1;
    if j > 3
        COP_avg(n,:) = mean(COP_avg(1:3,:));
        Room_avg(n,:) = mean(Room_avg(1:3,:));
        j = 1;
        n = n+1;
    end
end
hold off

COP_avg(6,:) = [];

subplot(3,1,2)
hold on
box on
for pp = 1:15
    plot(COP_avg(pp,:), 'r')
    plot(Room_avg(pp,:), 'r--')
end
hold off

subplot(3,1,3)
hold on
box on
COP_avg_total = mean(COP_avg);
Room_avg_total = mean(Room_avg);
plot(COP_avg_total,'r')
plot(Room_avg_total, 'r--')
hold off

COP_D075 = COP_avg;
COP_D075_total = COP_avg_total;
Room_D075 = Room_avg;
Room_D075_total = Room_avg_total;

clearvars -except COP_V075 COP_V075_total Room_V075 Room_V075_total COP_D075 COP_D075_total Room_D075 Room_D075_total

% Combination 0.75 Hz

cd('C:\Investigación\Moving Room\data')
load('ModalityCells.mat')

Fs = 50; % sampling frequency
Fc = 7.14; % cut-off frequency
order = 4;
j = 1; n = 1;

subplot(3,1,1)
figure(1)
hold on
box on
axis([0,1000, -1, 2])
for row = 1:48
    M = table2array(Combination_075_cell{row}(:,2:2:4));
    M_1 = [M(:,1),M(:,end)./1000];
    M_3 = M_1(M_1(:,2) > (5+1.34) & M_1(:,2) < (M_1(end,2)-5-1.34*1.5));
    COP = timeNormalize(M_3);

    % Get coefficients of the filter
    [b,a] = butter(order,Fc/(Fs/2),'low');
    COP_cf = filtfilt(b,a,COP);
    COP_n = normalize(COP_cf,'center');
    plot(COP_n, 'g')
    COP_avg(j,:) = COP_n;
    
    M_Room = [table2array(Combination_Room_075_cell{row}(:,4:4:8))];
    M_1_Room = [M_Room(:,1),M_Room(:,end)./1000];
    M_3_Room = M_1_Room(M_1_Room(:,2) > (5+1.34) & M_1_Room(:,2) < (M_1_Room(end,2)-5-1.34*1.5));
    Room = normalize(timeNormalize(M_3_Room));
    plot(Room, 'g--')
    Room_avg(j,:) = Room;
    
    clear M M_1 M_3 COP COP_cf COP_n
    j = j+1;
    if j > 3
        COP_avg(n,:) = mean(COP_avg(1:3,:));
        Room_avg(n,:) = mean(Room_avg(1:3,:));
        j = 1;
        n = n+1;
    end
end
hold off

COP_avg(6,:) = [];

subplot(3,1,2)
hold on
box on
for pp = 1:15
    plot(COP_avg(pp,:), 'g')
    plot(Room_avg(pp,:), 'g--')
end
hold off

subplot(3,1,3)
hold on
box on
COP_avg_total = mean(COP_avg);
Room_avg_total = mean(Room_avg);
plot(COP_avg_total,'g')
plot(Room_avg_total, 'g--')
hold off

COP_C075 = COP_avg;
COP_C075_total = COP_avg_total;
Room_C075 = Room_avg;
Room_C075_total = Room_avg_total;

Room_075_total = mean([Room_V075_total;Room_D075_total;Room_C075_total])./10;

clearvars -except COP_V075 COP_V075_total Room_V075 Room_V075_total COP_D075 COP_D075_total Room_D075 Room_D075_total COP_C075 COP_C075_total Room_C075 Room_C075_total Room_075_total
%
figure()
hold on
box on
title('COP 0.75 Hz')
plot(COP_V075_total, 'b', 'LineWidth', 1.5)
plot(COP_D075_total, 'r', 'LineWidth', 1.5)
plot(COP_C075_total, 'g', 'LineWidth', 1.5)
plot(Room_075_total, 'k', 'LineWidth', 1.5)
legend('Vision', 'Device', 'Combination','Wall', 'Location', 'SouthEast')
hold off

%%

clear all; close all; clc;
cd('C:\Investigación\Moving Room\data')
load('COP_avg.mat')

subplot(2,1,1)
hold on
box on
axis([0,1000, -.4, .4])
title('COP 0.10 Hz')
plot(COP_V01_total, 'b')
plot(COP_D01_total, 'r')
plot(COP_C01_total, 'g')
legend('Vision', 'Device', 'Combination', 'Location', 'SouthEast')
hold off

subplot(2,1,2)
hold on
box on
axis([0,1000, -.4, .4])
title('COP 0.75 Hz')
plot(COP_V075_total, 'b')
plot(COP_D075_total, 'r')
plot(COP_C075_total, 'g')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Average Plotting STILL (X)

%% Vision Still

clc; clear all; close all;
cd('C:\Investigación\Moving Room\data')
load('ModalityCells.mat')

Fs = 50; % sampling frequency
Fc = 7.14; % cut-off frequency
order = 4;
j = 1;

subplot(3,1,1)
figure(1)
hold on
box on
axis([0,1000, -1, 2])
for row = 1:16
    M = table2array(Vision_still_cell{row}(:,1:3:4));
    M_1 = [M(:,1),M(:,end)./1000];
    M_3 = M_1(M_1(:,2) > 5 & M_1(:,2) < 35);
    COPX = timeNormalize(M_3);

    % Get coefficients of the filter
    [b,a] = butter(order,Fc/(Fs/2),'low');
    COPX_cf = filtfilt(b,a,COPX);
    COPX_n = normalize(COPX_cf,'center');
    plot(COPX_n)

    COPX_avg(j,:) = COPX_n;

    clear M M_1 M_3 COP COP_cf COP_n
    j = j+1;
end
hold off

 COPX_avg(6,:) = [];

subplot(3,1,2)
hold on
box on
for pp = 1:15
    plot(COPX_avg(pp,:), 'k')
end
hold off

subplot(3,1,3)
hold on
box on
COPX_avg_total = mean(COPX_avg);
plot(COPX_avg_total)
hold off

COPX_VStill = COPX_avg;
COPX_VStill_total = COPX_avg_total;

clearvars -except COPX_VStill COPX_VStill_total
%
% Device 0.1 Hz

cd('C:\Investigación\Moving Room\data')
load('ModalityCells.mat')

Fs = 50; % sampling frequency
Fc = 7.14; % cut-off frequency
order = 4;
j = 1;

subplot(3,1,1)
figure(1)
hold on
box on
axis([0,1000, -1, 2])
for row = 1:16
    M = table2array(Device_still_cell{row}(:,1:3:4));
    M_1 = [M(:,1),M(:,end)./1000];
    M_3 = M_1(M_1(:,2) > 5 & M_1(:,2) < 35);
    COPX = timeNormalize(M_3);

    % Get coefficients of the filter
    [b,a] = butter(order,Fc/(Fs/2),'low');
    COPX_cf = filtfilt(b,a,COPX);
    COPX_n = normalize(COPX_cf,'center');
    plot(COPX_n)

    COPX_avg(j,:) = COPX_n;

    clear M M_1 M_3 COP COP_cf COP_n
    j = j+1;

end
hold off

COPX_avg(6,:) = [];

subplot(3,1,2)
hold on
box on
for pp = 1:15
    plot(COPX_avg(pp,:), 'k')
end
hold off

subplot(3,1,3)
hold on
box on
COPX_avg_total = mean(COPX_avg);
plot(COPX_avg_total)
hold off

COPX_DStill = COPX_avg;
COPX_DStill_total = COPX_avg_total;

clearvars -except COPX_VStill COPX_VStill_total COPX_DStill COPX_DStill_total

% Comb 0.1 Hz

cd('C:\Investigación\Moving Room\data')
load('ModalityCells.mat')

Fs = 50; % sampling frequency
Fc = 7.14; % cut-off frequency
order = 4;
j = 1;

subplot(3,1,1)
figure(1)
hold on
box on
axis([0,1000, -1, 2])
for row = 1:16
    M = table2array(Combination_still_cell{row}(:,1:3:4));
    M_1 = [M(:,1),M(:,end)./1000];
    M_3 = M_1(M_1(:,2) > 5 & M_1(:,2) < 35);
    COPX = timeNormalize(M_3);

    % Get coefficients of the filter
    [b,a] = butter(order,Fc/(Fs/2),'low');
    COPX_cf = filtfilt(b,a,COPX);
    COPX_n = normalize(COPX_cf,'center');
    plot(COPX_n)

    COPX_avg(j,:) = COPX_n;

    clear M M_1 M_3 COP COP_cf COP_n
    j = j+1;
end
hold off

COPX_avg(6,:) = [];

subplot(3,1,2)
hold on
box on
for pp = 1:15
    plot(COPX_avg(pp,:), 'k')
end
hold off

subplot(3,1,3)
hold on
box on
COPX_avg_total = mean(COPX_avg);
plot(COPX_avg_total)
hold off

COPX_CStill = COPX_avg;
COPX_CStill_total = COPX_avg_total;

clearvars -except COPX_VStill COPX_VStill_total COPX_DStill COPX_DStill_total COPX_CStill COPX_CStill_total

%%
figure()
hold on
box on
axis([0,1000, -.4, .4])
title('COP Still X')
plot(COPX_VStill_total, 'b')
plot(COPX_DStill_total, 'r')
plot(COPX_CStill_total, 'g')
legend('Vision', 'Device', 'Combination', 'Location', 'SouthEast')
hold off

%%

clear all; close all; clc;
cd('C:\Investigación\Moving Room\data')
load('Still.mat')

subplot(2,1,1)
hold on
box on
axis([0,1000, -.4, .4])
title('COP Y')
plot(COP_VStill_total, 'b')
plot(COP_DStill_total, 'r')
plot(COP_CStill_total, 'g')
legend('Vision', 'Device', 'Combination', 'Location', 'SouthEast')
hold off

subplot(2,1,2)
hold on
box on
axis([0,1000, -.4, .4])
title('COP X')
plot(COPX_VStill_total, 'b')
plot(COPX_DStill_total, 'r')
plot(COPX_CStill_total, 'g')
hold off
