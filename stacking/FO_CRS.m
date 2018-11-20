function [crs_section, crs_gathers, offsets_for_gathers, midpoints_for_gathers] ...
    = FO_CRS(traces, midpoints, offsets, dt, A, B, C, D, E, output_midpoints, output_offsets, aptx, apth)

%Output array
crs_section  = zeros(size(traces,1), length(output_midpoints), length(output_offsets));
n_traces = zeros(1, length(output_midpoints), length(output_offsets));
crs_gathers = {};
offsets_for_gathers = {};
midpoints_for_gathers = {};

t_vec = (1:size(traces,1))'.*dt;

xi =1;

LoopProgressReport('Making CRS stack')
for x0 = output_midpoints
    hi = 1;
    for h0 = output_offsets;

        a = collect_param_at_ind(A,xi,hi);
        b = collect_param_at_ind(B,xi,hi);
        c = collect_param_at_ind(C,xi,hi);
        d = collect_param_at_ind(D,xi,hi);
        e = collect_param_at_ind(E,xi,hi);
    
        %Select traces within aperture
        inds_in_apt = ...
            midpoints>= x0-aptx & ...
            midpoints<= x0+aptx & ...
            offsets>= h0-apth & ...
            offsets<= h0+apth;
        
        traces_in_apt = traces(:,inds_in_apt);
        midpoints_in_apt = midpoints(inds_in_apt);
        offsets_in_apt = offsets(inds_in_apt);

        gather = [];
        offset_for_gather = [];

        %loop through traces
        for ti = 1:size(traces_in_apt,2)
            h = offsets_in_apt(ti);
            x = offsets_in_apt(ti);

            trace = traces_in_apt(:,ti);

            %CRS moveout
            t_crs = sqrt( (t_vec + a.*(x-x0) + d.*(h-h0)).^2 + b.*(x-x0).^2 + c.*(h-h0).^2 + e.*(h-h0).*(x-x0) );

            %Linear interpolation
            time_corrected_trace = lininterp1f(t_vec, trace, t_crs, 0); %Linear interpolation

            %Add to output
            crs_section(:,xi,hi) = crs_section(:,xi,hi) + time_corrected_trace';
            n_traces(1,xi,hi) = n_traces(1,xi,hi) + 1;
            gather = [gather, time_corrected_trace'];
            offset_for_gather = [offset_for_gather, h0];
            midpoints_for_gathers = [midpoints_for_gathers, x0];

        end
        
    end
    crs_gathers{xi} = gather;
    offsets_for_gathers{xi} = offset_for_gather;
    xi = xi +1;
    LoopProgressReport(xi,length(output_midpoints))
end

%Normalize based on number of stacked traces
crs_section = crs_section./repmat(n_traces,size(crs_section,1),1);


end

function p = collect_param_at_ind(p,xi,hi);
    if length(size(p))==3;
        p = p(:, xi, hi);
    elseif length(size(p))==2 && size(p,2)~=1;
        p = p(:, xi);    
    else
        p = p(:);
    end
end
