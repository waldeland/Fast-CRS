function [crs_section,crs_gathers,offsets_for_gathers] = ZO_CRS(traces, midpoints, offsets, dt, A, B, C, output_midpoints, apt )
%Just a special case of FO CRS
output_offsets = [0];
apth = inf;
D = A*0;
E = A*0;
[crs_section,crs_gathers,offsets_for_gathers] = FO_CRS(traces, midpoints, offsets, dt, A, B, C, D, E, output_midpoints, output_offsets, apt, apth);
end


