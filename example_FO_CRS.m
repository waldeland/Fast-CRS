addpath('SegyMAT')
addpath('utils')
addpath('stacking')
addpath('semblance_search')
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

%% Full FO CRS parameter search + stack
%Select traces to output paramters for
%We only select one-offset to save time

output_offsets = [100];

%Define search ranges
param_search_range = struct();
param_search_range.A = linspace(-10^-3,10^-3,15);
param_search_range.B = linspace(-10^-6,10^-6,15);
param_search_range.C = linspace(-10^-7,10^-7,5);
param_search_range.D = linspace(-10^-4,10^-4,5);
param_search_range.E = linspace(-10^-6,10^-6,15);

% Search around velocity guide
[C,D] = convert_V(velocity_guide, output_offsets, dt);
initial_params = struct();
initial_params.A = 0;
initial_params.B = 0;
initial_params.C = C;
initial_params.D = D;
initial_params.E = 0;

%Define search aperture
aptx=100; apth=100; aptt=5;

%Do search
CRS_parameters = ...
    full_search_FO(traces, midpoints, offsets, dt,  output_midpoints, output_offsets, aptx, apth, aptt, param_search_range,  initial_params);

% Stack
[crs_section,crs_gathers,offsets_for_gathers] = FO_CRS(traces, midpoints, offsets, dt, initial_params, output_midpoints_, output_offsets_, apt, apth);
subplot(1,2,2);
imagesc(crs_gathers{end},imlim(crs_gathers{end}))
colormap('gray');
subplot(1,2,1);
imagesc(crs_section, imlim(crs_section))
colormap('gray');

%% Fast FO CRS parameter extraction + stack