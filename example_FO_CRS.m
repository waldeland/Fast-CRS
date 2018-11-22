addpath('SegyMAT')
addpath('utils')
addpath('utils/lininterp1f')

%% Load example data
data_path = 'example_data.segy';
[traces, H]=ReadSegy(data_path);

% Get parameters
offsets = [H.offset]./2;
dt = H(1).dt/1000/1000;
midpoints = [H.cdp].*25;
midpoints = midpoints-min(midpoints);%For simplicity we move the first x-location to 0

%Load velocity guide
velocity_guide = load('Vguide.mat');
velocity_guide = repmat(velocity_guide.Vguide ,1, length(unique(midpoints)));
output_midpoints = unique(midpoints);


%% NMO stack
[nmo_section,nmo_gathers,offsets_for_gathers] = NMO(traces, midpoints, offsets, dt, velocity_guide, output_midpoints);


subplot(1,2,2);
imagesc(nmo_gathers{end},imlim(nmo_gathers{end}))
colormap('gray');
subplot(1,2,1);
imagesc(nmo_section, imlim(nmo_section))
colormap('gray');


%% Full ZO CRS parameter search + stack
A =velocity_guide*0;
B = A;
C = 4./velocity_guide.^2;
apt = 100;

[crs_section,crs_gathers,offsets_for_gathers] = ZO_CRS(traces, midpoints, offsets, dt, A, B, C, output_midpoints, apt);
subplot(1,2,2);
imagesc(crs_gathers{end},imlim(crs_gathers{end}))
colormap('gray');
subplot(1,2,1);
imagesc(crs_section, imlim(crs_section))
colormap('gray');
%% Fast ZO CRS parameter extraction + stack

%% Full FO CRS parameter search + stack
%Select traces to output paramters for
%We only select some output-mid-points to save time
output_midpoints_ = output_midpoints(:,1:10);
velocity_guide_ = velocity_guide(:,1:10);
output_offsets_ = [0];

%Define search ranges
param_search_range = struct();
param_search_range.A = linspace(-10^-3,10^-3,5);
param_search_range.B = linspace(-10^-6,10^-6,5);
param_search_range.C_as_v = linspace(-200,200,5);
param_search_range.D_as_v = linspace(-200,200,5);
param_search_range.E = linspace(-10^-6,10^-6,5);

%
initial_params = struct();
initial_params.A = 0;
initial_params.B = 0;
initial_params.C = velocity_guide_;
initial_params.D = velocity_guide_;
initial_params.E = 0;

%Define search aperture
aptx=100; apth=100; aptt=5;

CRS_parameters = ...
    full_search_FO(traces, midpoints, offsets, dt,  output_midpoints_, output_offsets_, aptx, apth, aptt,param_search_range,  initial_params)

%% Fast FO CRS parameter extraction + stack