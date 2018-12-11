addpath('SegyMAT')
addpath('utils')
addpath('stacking')
addpath('semblance_search')
addpath('utils/lininterp1f')

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

%% NMO stack
[nmo_section,nmo_gathers,offsets_for_gathers] = NMO(traces, midpoints, offsets, dt, velocity_guide, output_midpoints);


subplot(1,2,2);
imagesc(nmo_gathers{100},imlim(nmo_gathers{100}))
colormap('gray');dx
subplot(1,2,1);
imagesc(nmo_section, imlim(nmo_section))
colormap('gray');


%% Full ZO CRS parameter search + stack

initial_params = struct();
initial_params.A = 0;
initial_params.B = 0;
initial_params.C = 4./velocity_guide.^2;

param_search_range = struct();
param_search_range.A = linspace(-10^-3,10^-3,15);
param_search_range.B = linspace(-10^-6,10^-6,15);
param_search_range.C = linspace(-10^-7,10^-7,11);

aptx=0;  aptt=5;

CRS_parameters = full_search_ZO(traces, midpoints, offsets, dt, output_midpoints, aptx, aptt, param_search_range,  initial_params);

apt = 100;
[crs_section,crs_gathers,offsets_for_gathers] = ZO_CRS(traces, midpoints, offsets, dt, CRS_parameters, output_midpoints, apt);

subplot(1,2,2);
imagesc(crs_gathers{end},imlim(crs_gathers{end}))
colormap('gray');
subplot(1,2,1);
imagesc(crs_section, imlim(crs_section))
colormap('gray');
%% Fast ZO CRS parameter extraction + stack
sigma_g = 1;
sigma_T = 5;
d = [dt,dx];
[A, B, coherency] = fastCRS(nmo_section, sigma_g, sigma_T, [dt,dx]);
    
CRS_parameters.A = A;
CRS_parameters.B = B;
CRS_parameters.C = 4./velocity_guide.^2;
apt = 25;
[crs_section,crs_gathers,offsets_for_gathers] = ZO_CRS(traces, midpoints, offsets, dt, CRS_parameters, output_midpoints, apt);

subplot(1,2,2);
imagesc(nmo_section,imlim(nmo_section))
colormap('gray');
subplot(1,2,1);
imagesc(crs_section, imlim(crs_section))
colormap('gray');

