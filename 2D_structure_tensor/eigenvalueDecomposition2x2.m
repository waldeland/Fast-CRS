function [e1,e2,v1,v2] = eigenvalueDecomposition2x2(T,sign_ind)

%Is sym
eps = .0000000001;
isSym = numel(T(:,:,1,1)).*eps > sum(sum(sum( abs( T(:,:,1,2) - T(:,:,2,1)) ))) ;
%Eigenvalue decomposition
if isSym
    [e1,e2,v1,v2] = eigRealSym2x2(T,1);
else
    
    % Output arrays
    out = T(:,:,1,1).*0;
    e1 = out; e2 = out;
    out = T(:,:,1:2).*0;
    v1 = out;   v2 = out;
    
    % Slow version
    for i = 1:size(T,1)
        for j = 1:size(T,2)
            
            t = reshape(T(i,j,:,:),2,2);
            [v,e] = eig(t);
            [e,ind] = sort(diag(e));
            e1(i,j) = e(2);
            e2(i,j) = e(1);
            %v = v(:,ind);
            
            if v(1,1) < 0
                v(:,1) = -v(:,1);
            end
            v1(i,j,:) = v(:,1);
            
            if v(1,2) < 0
                v(:,2) = -v(:,2);
            end
            v2(i,j,:) = v(:,2);
            
        end
    end
end
