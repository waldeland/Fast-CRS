

%Make toy example
[data, true_a, true_b, true_c, true_d, true_e, true_t, x, h, N] = toy_data_2D();
data = data+randn(size(data))*0.0001;

%% Proposed method
dt = 1; dx = 1; dy = 1;
[slopes, curvatures,coherency] = fastCRS(data,1,5, [dt,dx,dy]); 

%% make plot
%Function to get parameter values at surface
getP = @(param,ind) reshape( param(sub2ind(size(param),true_t(:),h(:)+N+1,x(:)+N+1,ones((N*2+1).^2,1).*ind )), 2*N+1,2*N+1);

%function to get color limits (symmetric around 0 and +- two times maximum)
getLim = @(param) [- max(abs(param(:))),  max(abs(param(:)))].*2;

subplot(3,5,1);
imagesc(true_a, getLim(true_a));
colorbar
axis equal tight off
title('True A')


subplot(3,5,2);
imagesc(true_d, getLim(true_d));
colorbar
axis equal tight off
title('True D')

subplot(3,5,3);
imagesc(true_b, getLim(true_b));
colorbar
axis equal tight off
title('True B')

subplot(3,5,4);
imagesc(true_c, getLim(true_c));
colorbar
axis equal tight off
title('True C')

subplot(3,5,5);
imagesc(true_e, getLim(true_e));
colorbar
axis equal tight off
title('True E')


subplot(3,5,6);
imagesc( getP(slopes,1), getLim(true_a));
colorbar
axis equal tight off
title('Estimated A')

subplot(3,5,7);
imagesc( getP(slopes,2), getLim(true_d));
colorbar
axis equal tight off
title('Estimated D')

subplot(3,5,8);
imagesc(getP(curvatures,1), getLim(true_b));
colorbar
axis equal tight off
title('Estimated B')

subplot(3,5,9);
imagesc(getP(curvatures,4), getLim(true_c));
colorbar
axis equal tight off
title('Estimated C')

subplot(3,5,10);
imagesc(getP(curvatures,2), getLim(true_e));
colorbar
axis equal tight off
title('Estimated E')

subplot(3,5,11);
imagesc(getP(coherency,1), [0,1]);
colorbar
axis equal tight off
title('Estimated Coherency')

colormap('jet')