function [mydepth] = getmydepth(fname,depth)

start

z0=-depth; % compute average above z0

nc=netcdf(fname);
pm=nc{'pm'}(:);
pn=nc{'pn'}(:);
h=nc{'h'}(:);
mask=nc{'mask_rho'}(:);
N=length(nc('s_rho'));
[M,L]=size(pm);
h = reshape(h,1,M,L);
h = repmat(h,[N 1 1]);
mask3d=reshape(mask,1,M,L);
mask3d=repmat(mask3d,[N 1 1]);   % Land = 0,  Sea = 1;
mask3d(mask3d == 0) = NaN;

ntime=length(nc('time'));
nstart=1;

for n = nstart:ntime
    zr=get_depths(fname,fname,n,'r');
    zr = zr.*mask3d;
    delta = abs(zr(2:end,:,:)-zr(1:end-1,:,:));
    delta = cat(1,h(1,:,:) - abs(zr(1,:,:)),delta);
    zr(zr>=z0) = 1;
    zr(zr<z0)=NaN;
    mydepth(n,:,:,:) = zr.*delta;
end

mydepth = permute(mydepth,[4,3,2,1]);







