function curveFull = interp1Segmented(curve,navIndices,linesPerShot,orient,method,navTiming,ACQTiming)

if nargin < 5
    method = 'pchip';
end
if nargin < 4
    orient = 'rows';
end
if ~strcmp(orient,'cols')
    curve = curve.';
end
if ndims(squeeze(curve)) == 1
    curve = curve(:);
end

navIndices = mod(navIndices-1,linesPerShot)+1;
recovery_starts = [1, find([0 diff(navIndices)]<0)];
curveFull = zeros(numel(recovery_starts)*linesPerShot,size(curve,2));

if nargin < 7
    navTiming = navIndices;
    ACQTiming = repmat(1:linesPerShot,[1 numel(recovery_starts)]);
end

for j = 1:numel(recovery_starts)-1 %each recovery curve separately
    t_ind  = recovery_starts(j):recovery_starts(j+1)-1; 
    Segidx  = navTiming(t_ind);
    Shotidx = ACQTiming((1:linesPerShot)+linesPerShot*(j-1));
    curveFull((j-1)*linesPerShot+(1:linesPerShot),:) = interp1(Segidx,curve(t_ind,:),Shotidx,method,'extrap');
end
if recovery_starts(end) < length(navIndices)
    t_ind  = recovery_starts(end):length(navIndices); 
    Segidx = navTiming(t_ind);
    Shotidx = ACQTiming((linesPerShot*j+1):end);
    curveFull((numel(recovery_starts)-1)*linesPerShot+(1:linesPerShot),:) = interp1(Segidx,curve(t_ind,:),Shotidx,method,'extrap');
end

if strcmp(orient,'rows')
    curveFull = curveFull.';
end
