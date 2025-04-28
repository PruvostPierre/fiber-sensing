path1 = 'birefringenceData1.mat';
path2 = 'birefringenceData2.mat';   
path3 = 'birefringenceData3.mat';
path4 = 'birefringenceData4.mat';
path5 = 'birefringenceData5.mat';
path6 = 'birefringenceGenData1.mat';
path7 = 'birefringenceGenData2.mat';
path8 = 'birefringenceGenData3.mat';
path9 = 'birefringenceGenData4.mat';
path10 = 'birefringenceGenData5.mat';

data_matrix1 = load(path1);
data1 = data_matrix1.birefringence; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
data1 = data1 - mean(data1, 3); % Subtract the mean along the third dimension

data_matrix2 =load(path2);
data2 = data_matrix2.birefringence; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
data2 = data2 - mean(data2, 3); % Subtract the mean along the third dimension

data_matrix3 =load(path3);
data3 = data_matrix3.birefringence; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
data3 = data3 - mean(data3, 3); % Subtract the mean along the third dimension

data_matrix4 =load(path4);
data4 = data_matrix4.birefringence; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
data4 = data4 - mean(data4, 3); % Subtract the mean along the third dimension

data_matrix5 =load(path5);
data5 = data_matrix5.birefringence; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
data5 = data5 - mean(data5, 3); % Subtract the mean along the third dimension

%combine all data above along third dimension
data = cat(3, data1, data2, data3, data4, data5); % Combine the matrices along the third dimension



dataGen_matrix1 =load(path6);
dataGen1 = dataGen_matrix1.birefringenceGen; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
dataGen1 = dataGen1 - mean(dataGen1, 3); % Subtract the mean along the third dimension

dataGen_matrix2 =load(path7);
dataGen2 = dataGen_matrix2.birefringenceGen; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
dataGen2 = dataGen2 - mean(dataGen2, 3); % Subtract the mean along the third dimension

dataGen_matrix3 =load(path8);
dataGen3 = dataGen_matrix3.birefringenceGen; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
dataGen3 = dataGen3 - mean(dataGen3, 3); % Subtract the mean along the third dimension

dataGen_matrix4 =load(path9);
dataGen4 = dataGen_matrix4.birefringenceGen; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
dataGen4 = dataGen4 - mean(dataGen4, 3); % Subtract the mean along the third dimension

dataGen_matrix5 =load(path10);
dataGen5 = dataGen_matrix5.birefringenceGen; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
dataGen5 = dataGen5 - mean(dataGen5, 3); % Subtract the mean along the third dimension

%combine all data above along third dimension
dataGen = cat(3, dataGen1, dataGen2, dataGen3, dataGen4, dataGen5); % Combine the matrices along the third dimension

%% Analysis of the data
%Compute and plot the standard deviation of the data along the third dimension
std_data = std(data, 0, 3); % Compute the standard deviation along the third dimension
std_dataGen = std(dataGen, 0, 3); % Compute the standard deviation along the third dimension



% Custom red-white-blue colormap
n = 256;
rwb = [linspace(0,1,n/2)', linspace(0,1,n/2)', ones(n/2,1); ...
       ones(n/2,1), linspace(1,0,n/2)', linspace(1,0,n/2)'];

figure;

% First subplot: Applied Birefringence Heat Map
subplot(1, 2, 1);
disp('Birefringence heat map');
x_axis_meters = (1:2927) * 1.025; % Convert reflector indices to meters
y_axis_seconds = (1:495) * 2.62*1e-3; % Convert code indices to seconds
birefringenceHeatMap = squeeze(dataGen(:,:,1:end)); % Extract birefringence data
imagesc(x_axis_meters, y_axis_seconds, abs(birefringenceHeatMap.')); % Transpose for correct orientation
colorbar;
%colormap(rwb); % Apply custom colormap
xlabel('Distance (m)', 'FontSize', 14);
ylabel('Time (s)', 'FontSize', 14);
xlim([1300 1850]);
%title('Applied Birefringence Heat Map');
grid on;

% Second subplot: Estimated Birefringence Heat Map
subplot(1, 2, 2);
disp('Birefringence heat map');
x_axis_meters = (1:292) * 10.25; % Convert reflector indices to meters
y_axis_seconds = (1:495) * 2.62*1e-3; % Convert code indices to seconds
birefringenceHeatMap = squeeze(data(:,:,1:end)); % Extract birefringence data
imagesc(x_axis_meters, y_axis_seconds, abs(birefringenceHeatMap.')); % Transpose for correct orientation
colorbar;
%colormap(rwb); % Apply custom colormap
xlabel('Distance (m)', 'FontSize', 14);
ylabel('Time (s)', 'FontSize', 14);
xlim([1300, 1850]);
%title('Estimated Birefringence Heat Map');
grid on;



figure;
x_axis_meters = (1:292) * 10.25; % Convert reflector indices to meters
stdBirefringenceSelect = std(abs(data(1, :, 1:end-1)), 0, 3); % Compute std across time (3rd dimension), discarding the first 20 codes
% Create a subplot for the two standard deviation plots
subplot(1, 2, 2);
plot(x_axis_meters, stdBirefringenceSelect, '-o');
xlabel('Distance (m)', 'FontSize', 14);
ylabel('Standard deviation of Estimated Birefringence (radians)', 'FontSize', 12);
%title('Standard Deviation of Estimated Birefringence Across Time');
grid on;
subplot(1, 2, 1);
x_axis_meters = (1:2926) * 1.025; % Convert reflector indices to meters
stdBirefringenceGen = std(abs(dataGen(1, :, 1:end-1)), 0, 3); % Compute std across time (3rd dimension), discarding the first 20 codes*
disp(size(dataGen));
plot(x_axis_meters, stdBirefringenceGen, '-o', 'Color', '#A2142F');
xlabel('Distance (m)', 'FontSize', 14);
ylabel('Standard Deviation of Applied Birefringence (radians)', 'FontSize', 12);
%title('Standard Deviation of Applied Birefringence Across Time');
grid on;


%%visualization of generated data
% First subplot: Applied Birefringence Heat Map
figure;
subplot(1, 2, 1);
disp('Birefringence heat map');
x_axis_meters = (1:2927) * 1.025; % Convert reflector indices to meters
y_axis_seconds = (1:495) * 2.62*1e-3; % Convert code indices to seconds
birefringenceHeatMap = squeeze((dataGen(:,1463:1610,1:end))); % Extract birefringence data
imagesc(1500:1650, y_axis_seconds, birefringenceHeatMap.'); % Transpose for correct orientation
colorbar;
colormap(rwb); % Apply custom colormap
xlabel('Distance (m)', 'FontSize', 14);
ylabel('Time (s)', 'FontSize', 14);
%xlim([1300 1850]);
%title('Applied Birefringence Heat Map');
grid on;

disp('Birefringence heat map');
subplot(1, 2, 2);
x_axis_meters = (1:292) * 10.25; % Convert reflector indices to meters
y_axis_seconds = (1:495) * 2.62*1e-3; % Convert code indices to seconds
birefringenceHeatMap = squeeze((data(:,146:161,1:end))); % Extract birefringence data
imagesc(1500:1650, y_axis_seconds, birefringenceHeatMap.'); % Transpose for correct orientation
colorbar;
colormap(rwb); % Apply custom colormap
xlabel('Distance (m)', 'FontSize', 14);
ylabel('Time (s)', 'FontSize', 14);
%xlim([1300 1850]);
%title('Applied Birefringence Heat Map');
grid on;


%% 2304 data
path1 = 'birefringenceData2304_1.mat';
path2 = 'birefringenceData2304_2.mat';   
path3 = 'birefringenceData2304_3.mat';
path4 = 'birefringenceData2304_4.mat';
path5 = 'birefringenceData2304_5.mat';
path12 = 'birefringenceData2304_6.mat';
path6 = 'birefringenceGenData2304_1.mat';
path7 = 'birefringenceGenData2304_2.mat';
path8 = 'birefringenceGenData2304_3.mat';
path9 = 'birefringenceGenData2304_4.mat';
path10 = 'birefringenceGenData2304_5.mat';
path11 = 'birefringenceGenData2304_6.mat';

data_matrix1 = load(path1);
data1 = data_matrix1.birefringence; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
data1 = data1 - mean(data1, 3); % Subtract the mean along the third dimension

data_matrix2 =load(path2);
data2 = data_matrix2.birefringence; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
data2 = data2 - mean(data2, 3); % Subtract the mean along the third dimension

data_matrix3 =load(path3);
data3 = data_matrix3.birefringence; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
data3 = data3 - mean(data3, 3); % Subtract the mean along the third dimension

data_matrix4 =load(path4);
data4 = data_matrix4.birefringence; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
data4 = data4 - mean(data4, 3); % Subtract the mean along the third dimension

data_matrix5 =load(path5);
data5 = data_matrix5.birefringence; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
data5 = data5 - mean(data5, 3); % Subtract the mean along the third dimension

data_matrix6 =load(path12);
data6 = data_matrix6.birefringence; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
data6 = data6 - mean(data6, 3); % Subtract the mean along the third dimension


%combine all data above along third dimension
data = cat(3, data1, data2, data3, data4, data5, data6); % Combine the matrices along the third dimension


dataGen_matrix1 =load(path6);
dataGen1 = dataGen_matrix1.birefringenceGen; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
dataGen1 = dataGen1 - mean(dataGen1, 3); % Subtract the mean along the third dimension

dataGen_matrix2 =load(path7);
dataGen2 = dataGen_matrix2.birefringenceGen; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
dataGen2 = dataGen2 - mean(dataGen2, 3); % Subtract the mean along the third dimension

dataGen_matrix3 =load(path8);
dataGen3 = dataGen_matrix3.birefringenceGen; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
dataGen3 = dataGen3 - mean(dataGen3, 3); % Subtract the mean along the third dimension

dataGen_matrix4 =load(path9);
dataGen4 = dataGen_matrix4.birefringenceGen; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
dataGen4 = dataGen4 - mean(dataGen4, 3); % Subtract the mean along the third dimension

dataGen_matrix5 =load(path10);
dataGen5 = dataGen_matrix5.birefringenceGen; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
dataGen5 = dataGen5 - mean(dataGen5, 3); % Subtract the mean along the third dimension

dataGen_matrix6 =load(path11);
dataGen6 = dataGen_matrix6.birefringenceGen; % Extract the 'data' field from the loaded structure
%take data-mean(data) along third dimension
dataGen6 = dataGen6 - mean(dataGen6, 3); % Subtract the mean along the third dimension

%combine all data above along third dimension
dataGen = cat(3, dataGen1, dataGen2, dataGen3, dataGen4, dataGen5, dataGen6); % Combine the matrices along the third dimension

%% Analysis of the data
%Compute and plot the standard deviation of the data along the third dimension
std_data = std(data, 0, 3); % Compute the standard deviation along the third dimension
std_dataGen = std(dataGen, 0, 3); % Compute the standard deviation along the third dimension



% Custom red-white-blue colormap
n = 256;
rwb = [linspace(0,1,n/2)', linspace(0,1,n/2)', ones(n/2,1); ...
       ones(n/2,1), linspace(1,0,n/2)', linspace(1,0,n/2)'];

figure;

% First subplot: Applied Birefringence Heat Map
subplot(1, 2, 1);
disp('Birefringence heat map');
x_axis_meters = (1:2927) * 1.025; % Convert reflector indices to meters
y_axis_seconds = (1:495) * 2.62*1e-3; % Convert code indices to seconds
birefringenceHeatMap = squeeze(dataGen(:,:,1:end)); % Extract birefringence data
imagesc(x_axis_meters, y_axis_seconds, abs(birefringenceHeatMap.')); % Transpose for correct orientation
colorbar;
%colormap(rwb); % Apply custom colormap
xlabel('Distance (m)', 'FontSize', 14);
ylabel('Time (s)', 'FontSize', 14);
xlim([1300 1850]);
%title('Applied Birefringence Heat Map');
grid on;

% Second subplot: Estimated Birefringence Heat Map
subplot(1, 2, 2);
disp('Birefringence heat map');
x_axis_meters = (1:292) * 10.25; % Convert reflector indices to meters
y_axis_seconds = (1:495) * 2.62*1e-3; % Convert code indices to seconds
birefringenceHeatMap = squeeze(data(:,:,1:end)); % Extract birefringence data
imagesc(x_axis_meters, y_axis_seconds, abs(birefringenceHeatMap.')); % Transpose for correct orientation
colorbar;
%colormap(rwb); % Apply custom colormap
xlabel('Distance (m)', 'FontSize', 14);
ylabel('Time (s)', 'FontSize', 14);
xlim([1300, 1850]);
%title('Estimated Birefringence Heat Map');
grid on;



figure;
x_axis_meters = (1:292) * 10.25; % Convert reflector indices to meters
stdBirefringenceSelect = std(abs(data(1, :, 1:end-1)), 0, 3); % Compute std across time (3rd dimension), discarding the first 20 codes
% Create a subplot for the two standard deviation plots
subplot(1, 2, 2);
plot(x_axis_meters, stdBirefringenceSelect, '-o');
xlabel('Distance (m)', 'FontSize', 14);
ylabel('Standard deviation of Estimated Birefringence (radians)', 'FontSize', 12);
%title('Standard Deviation of Estimated Birefringence Across Time');
grid on;
subplot(1, 2, 1);
x_axis_meters = (1:2926) * 1.025; % Convert reflector indices to meters
stdBirefringenceGen = std(abs(dataGen(1, :, 1:end-1)), 0, 3); % Compute std across time (3rd dimension), discarding the first 20 codes*
disp(size(dataGen));
plot(x_axis_meters, stdBirefringenceGen, '-o', 'Color', '#A2142F');
xlabel('Distance (m)', 'FontSize', 14);
ylabel('Standard Deviation of Applied Birefringence (radians)', 'FontSize', 12);
%title('Standard Deviation of Applied Birefringence Across Time');
grid on;


%%visualization of generated data
% First subplot: Applied Birefringence Heat Map
figure;
subplot(1, 2, 1);
disp('Birefringence heat map');
x_axis_meters = (1:2927) * 1.025; % Convert reflector indices to meters
y_axis_seconds = (1:495) * 2.62*1e-3; % Convert code indices to seconds
birefringenceHeatMap = squeeze((dataGen(:,1463:1610,1:end))); % Extract birefringence data
imagesc(1500:1650, y_axis_seconds, birefringenceHeatMap.'); % Transpose for correct orientation
colorbar;
colormap(rwb); % Apply custom colormap
xlabel('Distance (m)', 'FontSize', 14);
ylabel('Time (s)', 'FontSize', 14);
%xlim([1300 1850]);
%title('Applied Birefringence Heat Map');
grid on;

disp('Birefringence heat map');
subplot(1, 2, 2);
x_axis_meters = (1:292) * 10.25; % Convert reflector indices to meters
y_axis_seconds = (1:495) * 2.62*1e-3; % Convert code indices to seconds
birefringenceHeatMap = squeeze((data(:,146:161,1:end))); % Extract birefringence data
imagesc(1500:1650, y_axis_seconds, birefringenceHeatMap.'); % Transpose for correct orientation
colorbar;
colormap(rwb); % Apply custom colormap
xlabel('Distance (m)', 'FontSize', 14);
ylabel('Time (s)', 'FontSize', 14);
%xlim([1300 1850]);
%title('Applied Birefringence Heat Map');
grid on;