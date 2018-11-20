%% stackNMO
%==========================================================================
% Stack data according to NMO correction
%==========================================================================
function [nmo_section,nmo_gathers,offsets_for_gathers] = NMO(traces, midpoints, offsets, dt, velocity_guide, output_midpoints)

%Output array
nmo_section  = zeros(size(traces,1), length(output_midpoints));
n_traces = zeros(1, length(output_midpoints));
nmo_gathers = {};
offsets_for_gathers = {};

t_vec = (1:size(traces,1))'.*dt;

xi =1;
for x = output_midpoints
    v = velocity_guide(:,xi);
    
    %Select traces with given midpoint
    inds_for_midpoint = midpoints==x;
    traces_for_midpoint = traces(:,inds_for_midpoint);
    offsets_for_traces = offsets(inds_for_midpoint);
    
    gather = [];
    offset_for_gather = [];
    
    %loop through traces
    for ti = 1:size(traces_for_midpoint,2)
        h = offsets_for_traces(ti);
        
        trace = traces_for_midpoint(:,ti);
        
        %NMO moveout
        t_nmo = sqrt( t_vec.^2 + 4*h.^2./v.^2);
        
        %Linear interpolation
        time_corrected_trace = lininterp1f(t_vec, trace, t_nmo, 0); %Linear interpolation
        
        %Add to output
        nmo_section(:,xi) = nmo_section(:,xi) + time_corrected_trace';
        n_traces(1,xi) = n_traces(1,xi) + 1;
        gather = [gather, time_corrected_trace'];
        offset_for_gather = [offset_for_gather, h];
    
    end
    
    nmo_gathers{xi} = gather;
    offsets_for_gathers{xi} = offset_for_gather;
    xi = xi +1;
end

%Normalize based on number of stacked traces
pnmo_section = nmo_section./repmat(n_traces,size(nmo_section,1),1);


end

