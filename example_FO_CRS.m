addpath('SegyMAT')
addpath('utils')
addpath('stacking')
addpath('semblance_search')
addpath('utils/lininterp1f')
addpath('nD_structure_tensor')

%% Load example data
data_path = 'example_data.segy'; 
[traces, H]=ReadSegy(data_path);
dx = 12.5;

% Get/set parameters
offsets = [H.offset]./2;
dt = H(1).dt/1000/1000;
midpoints = [H.cdp].*dx;
output_midpoints = unique(midpoints);


%Load velocity guide
velocity_guide = load('velocity_guide.mat');
velocity_guide = velocity_guide.velocity_guide;

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
[crs_section,crs_gathers,offsets_for_gathers] = FO_CRS(traces, midpoints, offsets, dt, initial_params, output_midpoints, output_offsets, aptx, apth);
subplot(1,2,2);
imagesc(crs_gathers{end},imlim(crs_gathers{end}))
colormap('gray');
subplot(1,2,1);
imagesc(crs_section, imlim(crs_section))
colormap('gray');

%% Fast FO CRS parameter extraction + stack
output_offsets = [100];

%Step 1: Convert V to C for the given offset
[C,D] = convert_V(velocity_guide, output_offsets, dt);


%Step 2: Create a mini-stack with aptx=0
initial_params = struct();
initial_params.A = 0;
initial_params.B = 0;
initial_params.C = C;
initial_params.D = D;
initial_params.E = 0;
aptx=0;
apth=50;
[crs_section,crs_gathers,offsets_for_gathers] = FO_CRS(traces, midpoints, offsets, dt, initial_params, output_midpoints, output_offsets, aptx, apth);

%Step 3: Estimate A and B
sigma_g = 1;
sigma_T = 5;
d = [dt,dx];
[A, B, coherency] = fastCRS(crs_section, sigma_g, sigma_T, [dt,dx]);

%Step 4: Estimate E with semblance search
aptx=50;
apth=50;
aptt=5;

initial_params.A = A;
initial_params.B = B;
initial_params.E = 0;

param_search_range = struct();
param_search_range.A = [0];
param_search_range.B = [0];
param_search_range.C = [0];
param_search_range.D = [0];
param_search_range.E = linspace(-10^-6,10^-6,15);

CRS_parameters = ...
    full_search_FO(traces, midpoints, offsets, dt,  output_midpoints, output_offsets, aptx, apth, aptt, param_search_range,  initial_params);
    
% Stack
[crs_section,crs_gathers,offsets_for_gathers] = FO_CRS(traces, midpoints, offsets, dt, initial_params, output_midpoints, output_offsets, aptx, apth);
subplot(1,2,2);
imagesc(crs_gathers{end},imlim(crs_gathers{end}))
colormap('gray');
subplot(1,2,1);
imagesc(crs_section, imlim(crs_section))
colormap('gray');