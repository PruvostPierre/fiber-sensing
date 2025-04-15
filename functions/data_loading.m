function data=data_loading(path)
% Load data from the specified path
data_matrix =load(path, 'data');
data = data_matrix.data; % Extract the 'data' field from the loaded structure
%time_sampling = 400;
%gauge = 3.2;
%x = (0:size(data, 1)-1) * gauge * 1e-3; %*gauge
%t = (0:size(data, 2)-1) / time_sampling;
