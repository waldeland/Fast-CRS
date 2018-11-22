function CRS_parameters = ...
    full_search_FO(traces, midpoints, offsets, dt, ...
    output_midpoints, output_offsets, aptx, apth, aptt, param_search_range,  initial_params)

%Array to keep track of the best parameter-combination
highest_semblance = zeros(size(traces,1), length(output_midpoints), length(output_offsets));

%Output arrays
best_a = zeros(size(traces,1), length(output_midpoints), length(output_offsets));
best_b = zeros(size(traces,1), length(output_midpoints), length(output_offsets));
best_c = zeros(size(traces,1), length(output_midpoints), length(output_offsets));
best_d = zeros(size(traces,1), length(output_midpoints), length(output_offsets));
best_e = zeros(size(traces,1), length(output_midpoints), length(output_offsets));

%Loop through parameter combinations
LoopProgressReport('Running FO CRS search'); i = 1; 
N = length(param_search_range.B)*length(param_search_range.A)*length(param_search_range.C_as_v)*length(param_search_range.D_as_v)*length(param_search_range.E);

for A = param_search_range.A
    A = A + initial_params.A;
    
    for B = param_search_range.B
        B = B + initial_params.B;
        
        for C = param_search_range.C_as_v
            c_as_v = sqrt(4./initial_params.C) + C;
            C = 4./c_as_v.^2;
            
            for D = param_search_range.D_as_v
                d_as_v = sqrt(4./initial_params.D) + D;
                D = 4./d_as_v.^2;
                
                for E = param_search_range.E
                    E = E + initial_params.E;
                    
                    %Loop through x,h points in output-sections
                    xi=0;
                    for x0 = output_midpoints
                        xi = xi+1;
                        hi = 0;
                        for h0 = output_offsets;
                            hi = hi+1;
                            
                            %Select traces inside aperture
                            inds_in_apt = ...
                                midpoints>= x0-aptx & ...
                                midpoints<= x0+aptx & ...
                                offsets>= h0-apth & ...
                                offsets<= h0+apth;
                            
                            %Loop through time samples
                            for ti = 1:size(traces,1);
                                t0 = ti*dt;
                                
                                %Select parameters to test for index
                                a = collect_param_at_ind(A,ti,xi,hi);
                                b = collect_param_at_ind(B,ti,xi,hi);
                                c = collect_param_at_ind(C,ti,xi,hi);
                                d = collect_param_at_ind(D,ti,xi,hi);
                                e = collect_param_at_ind(E,ti,xi,hi);
                                
                                traces_in_apt = traces(:,inds_in_apt);
                                midpoints_in_apt = midpoints(inds_in_apt);
                                offsets_in_apt = offsets(inds_in_apt);
                                
                                [semb] = get_semblance_at_location_for_parameters(...
                                    traces_in_apt, a, b, c, d, e, midpoints_in_apt, offsets_in_apt, ...
                                    h0, x0, t0, aptt, dt);
                                
                                %Store parameters if they gave better
                                %results
                                if semb>highest_semblance(ti,xi,hi)
                                    highest_semblance(ti,xi,hi) = semb;
                                    best_a(ti,xi,hi) = a;
                                    best_b(ti,xi,hi) = b;
                                    best_c(ti,xi,hi) = c;
                                    best_d(ti,xi,hi) = d;
                                    best_e(ti,xi,hi) = e;
                                end
                            end
                        end
                    end
                    LoopProgressReport(i,N);i = i+1;
                end
            end
        end
    end 
end
%Put in struct
CRS_parameters = struct();
CRS_parameters.A = best_a;
CRS_parameters.B = best_b;
CRS_parameters.C = best_c;
CRS_parameters.D = best_d;
CRS_parameters.E = best_e;
end

function p = collect_param_at_ind(p,ti,xi,hi);
    if length(size(p))==3;
        p = p(ti, xi, hi);
    elseif length(size(p))==2 && size(p,2)~=1;
        p = p(ti, xi);    
    elseif length(size(p))==2 && size(p,2)==1 && size(p,1)~=1;
        p = p(ti);    
    elseif length(size(p))==1;
        p = p(ti);    
    else
        p = p;
    end
end



function [semb] = get_semblance_at_location_for_parameters(...
    traces_in_apt, a, b, c, d, e, ...
    midpoints_in_apt, offsets_in_apt, h0, x0, t0, aptt, dt)

Napt = length(midpoints_in_apt);
Nt = size(traces_in_apt,1);
%Corrected amplitues
amp = zeros(aptt*2+1, Napt);

%Loop through traces inside aperture in midpoints/offsets
N = 0;
for ia = 1:Napt;
    h = offsets_in_apt(ia);
    x = midpoints_in_apt(ia);
    
    %Put part of the CRS-formula outside the loop to save some computation
    %overhead
    
    dm = x-x0;
    a_dm = a.*dm;
    b_dm_2 = b.*dm.^2;
    
    dh = h-h0;
    ch2 = c.*dh.*dh ;
    d_dh = d.*dh;
    
    ehm = e.*dh.*dm ;
    
    %Aperture in t
    for delta_it = -aptt:aptt
        
        t = t0 + delta_it*dt;
        
        t_crs = sqrt( ( t  + a_dm +  d_dh).^2 + b_dm_2 + ch2 + ehm );
        
        %Linear interpolation
        tmp = (t_crs-t0)/dt;
        upp = ceil( tmp );
        low = floor( tmp );
        
        %Check bounds
        if upp <= Nt & low > 1
            weigth = ( tmp - low);
            amp(delta_it+aptt+1, ia) = traces_in_apt(upp,ia)*weigth + traces_in_apt(low,ia)*(1-weigth);
            N = N + 1;
        end
    end
end
%Compute semblance
semb = sum( sum(amp, 2).^2 ,1) ./ sum(amp(:).^2)./N;
end