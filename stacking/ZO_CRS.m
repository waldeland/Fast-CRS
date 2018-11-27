function [crs_section,crs_gathers,offsets_for_gathers] = ZO_CRS(traces, midpoints, offsets, dt, CRS_parameters, output_midpoints, apt )
%Just a special case of FO CRS
output_offsets = [0];
apth = inf;
CRS_parameters.D = CRS_parameters.A*0;
CRS_parameters.E = CRS_parameters.A*0;
[crs_section,crs_gathers,offsets_for_gathers] = FO_CRS(traces, midpoints, offsets, dt, CRS_parameters, output_midpoints, output_offsets, apt, apth);
end


