function [C,D] = convert_V(V,offsets, dt)
%Convert V to C0 (c at offset=0)
C_Z0 = 4./V.^2;

C = zeros(size(V,1),size(V,2),length(offsets));
D = zeros(size(V,1),size(V,2),length(offsets));

t_0 = (1:size(V,1))'*dt; %Time at a central point
t_0 = repmat(t_0,1,size(V,2));


i=0;
for offset = offsets
    i = i+1;
    
    % Get t_Z0 and C_Z0 for the given offset
    t_Z0 = invNMO(t_0, V, offset, dt );
    C_Z0 = invNMO(C_Z0, V, offset, dt );

    %When there is no point in ZO for a given central point
    region_outside = t_Z0 == 0;
    t_Z0(region_outside) = min(t_Z0(~region_outside));
    C_Z0(region_outside) = max(C_Z0(~region_outside));

    % Comute D and C
    D(:,:,i) = C_Z0.*offset./t_0;
    C(:,:,i) = t_Z0.^2.*C_Z0./(t_0.^2);
end
end


%Inverse NMO shift
function out_data = invNMO(data,V,offset,dt)
out_data  = zeros(size(data));


t_0 = (1:size(V,1))'*dt;

%Loop trhough traces and do inverse NMO
for ix = 1:size(V,2)
    trace = data(:,ix);
        
    t_nmo = real(sqrt( t_0.^2 - 4*offset.^2./V(:,ix).^2));
    time_corrected_trace = lininterp1f(t_0 , trace, t_nmo, 0); %Linear interpolation with nan (4th argument) outside valid range
   
    out_data(:,ix) = time_corrected_trace;
end

out_data(isnan(out_data))=0;
end