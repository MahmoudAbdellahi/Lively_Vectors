function squeezed = lv_squeeze(z, squeeze_dims) % z is mdim tensor and dims to squeeze
% squeezes specific dimensions in the data determined by the vector squeeze_dims..
% example:
% x = [1 2 3 ; 4 5 6];
% y(:,:,2) =  x ;
% z(1,1,:,:,1) = x;
 
sz = size(z); if any(squeeze_dims>1), error('lv: you are squeezing a dimension that is non-singleton which will corrupt data!!'); end
sz(squeeze_dims)=[];
squeezed = reshape(z, sz);


end