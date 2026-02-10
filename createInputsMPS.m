function [di, partialTi, ki, sp] = createInputsMPS(lulcDir,mps,maskDir,sizeMap)

% CREATE MASK
% Mask with shape of image
if mps.useMask == true
    mask = single(readgeoraster(maskDir));
else
    mask = ones(sizeMap);
end
sg = nan(size(mask));

% CREATE LATLON LAYERS
% Generate longitude (lon) and latitude (lat) maps
lon = zeros(size(sg));
lat = zeros(size(sg));
% Not real latlon, just info of where in image
lon_start = 0;  % Starting longitude
lat_start = 0;  % Starting latitude
% Populate lon and lat based on pixel positions
for la = 1:size(sg,1)
    for lo = 1:size(sg,2)
        lon(la,lo) = lon_start + lo;
        lat(la,lo) = lat_start + la;
    end
end
lat(mask==0) = nan;
lon(mask==0) = nan;
latlon = cat(3, lat, lon);

% LOAD LULC MAP
lulc = single(readgeoraster(lulcDir));
lulc(mask==0) = nan;

% CREATE PARTIAL TI
partialTi = cat(3, latlon, lulc);

% CREATE DI
di = cat(3, sg, partialTi);

% CREATE KI
ki = ones(mps.kernel_dims);
for w = 1:length(mps.kernel_weights)
    ki(:,:,w) = ki(:,:,w) .* mps.kernel_weights(w);
end

% CREATE SP
if mps.useSP == true
    rng(mps.seed)
    values = randperm(numel(mask));
    spSingle = reshape(values, size(mask));
    if mps.useMask == true
        filler = ones(size(partialTi)) * -inf;
        spSingle(mask==0) = -inf;
        sp = cat(3,spSingle,filler);
    end
else
    if mps.useMask == true
        di(isnan(di)) = -999;
        di(mask==1)   = nan;
    end
    sp = [];
end

end
