function CRS_parameters = ...
    full_search_ZO(traces, midpoints, offsets, dt, ...
    output_midpoints, aptx, aptt, param_search_range,  initial_params)

%Array to keep track of the best parameter-combination
highest_semblance = zeros(size(traces,1), length(output_midpoints));

%Output arrays
best_a = zeros(size(traces,1), length(output_midpoints));
best_b = zeros(size(traces,1), length(output_midpoints));
best_c = zeros(size(traces,1), length(output_midpoints));


%Loop through parameter combinations
LoopProgressReport('Running FO CRS search'); i = 1;
N = length(param_search_range.B)*length(param_search_range.A)*length(param_search_range.C);

for A = param_search_range.A
    A = A + initial_params.A;
    
    for B = param_search_range.B
        B = B + initial_params.B;
        
        for C = param_search_range.C
            C = C + initial_params.C;
            
            
            %Loop through x,h points in output-sections
            xi=0;
            for x0 = output_midpoints
                xi = xi+1;
                
                
                %Select traces inside aperture
                inds_in_apt = ...
                    midpoints>= x0-aptx & ...
                    midpoints<= x0+aptx;
                
                %Loop through time samples
                for ti = 1:size(traces,1);
                    t0 = ti*dt;
                    
                    %Select parameters to test for index
                    a = collect_param_at_ind(A,ti,xi);
                    b = collect_param_at_ind(B,ti,xi);
                    c = collect_param_at_ind(C,ti,xi);
                    
                    traces_in_apt = traces(:,inds_in_apt);
                    midpoints_in_apt = midpoints(inds_in_apt);
                    offsets_in_apt = offsets(inds_in_apt);
                    
                    [semb] = get_semblance_at_location_for_parameters(...
                        traces_in_apt, a, b, c, midpoints_in_apt,offsets_in_apt, ...
                        x0, t0, aptt, dt);
                    
                    %Store parameters if they gave better
                    %results
                    if semb>highest_semblance(ti,xi)
                        highest_semblance(ti,xi) = semb;
                        best_a(ti,xi) = a;
                        best_b(ti,xi) = b;
                        best_c(ti,xi) = c;
                        
                    end
                    

                end
            end
            LoopProgressReport(i,N);i = i+1;
        end
    end
end
%Put in struct
CRS_parameters = struct();
CRS_parameters.A = best_a;
CRS_parameters.B = best_b;
CRS_parameters.C = best_c;

end

function p = collect_param_at_ind(p,ti,xi)
if length(size(p))==3;
    p = p(ti, xi);
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
    traces_in_apt, a, b, c,  ...
    midpoints_in_apt, offsets_in_apt,  x0, t0, aptt, dt)

Napt = length(midpoints_in_apt);
Nt = size(traces_in_apt,1);
%Corrected amplitues
amp = zeros(aptt*2+1, Napt);

%Loop through traces inside aperture in midpoints/offsets
N = 0;
for ia = 1:Napt;
    
    x = midpoints_in_apt(ia);
    h = offsets_in_apt(ia);
    
    %Put part of the CRS-formula outside the loop to save some computation
    %overhead
    
    dm = x-x0;
    a_dm = a.*dm;
    b_dm_2 = b.*dm.^2;
    
    dh = h;
    ch2 = c.*dh.*dh ;

    %Aperture in t
    for delta_it = -aptt:aptt
        
        t = t0 + delta_it*dt;
        
        t_crs = sqrt( ( t  + a_dm ).^2 + b_dm_2 + ch2 );
        
        %Linear interpolation
        tmp = (t_crs)/dt;
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
