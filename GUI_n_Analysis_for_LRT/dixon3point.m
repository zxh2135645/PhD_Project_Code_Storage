function [W,F,fB0,R2s] = dixon3point(S, TEs, varargin)
% DIXON3POINT  Three-point Dixon water–fat separation (single-peak fat).
%
% Inputs
%   S    : [Ny Nx 3] or [Ny Nx Necho] complex images (Necho>=3; uses first 3)
%   TEs  : [1xNecho] echo times in seconds (use first 3)
%
% Name-Value (all optional)
%   'DeltaFHz'   : fat–water chemical shift in Hz (default 428 for ~3T)
%   'RefineB0'   : logical, do 1-D Gauss-Newton refinement (default true)
%   'DoR2s'      : logical, estimate R2* (default false, more stable off)
%   'B0RangeHz'  : scalar, clip initial fB0 to +/- range (default 300)
%   'NumIter'    : integer, B0 refinement iterations (default 3)
%   'MedianWin'  : odd int for median filter on fB0 (default 3; set 0 to skip)
%
% Outputs
%   W,F  : complex water and fat images
%   fB0  : B0 field map [Hz]
%   R2s  : R2* [1/s] (zeros if DoR2s=false)
%
% Example:
%   [W,F,fB0] = dixon3point(S, [1.2e-3 2.4e-3 3.6e-3], 'DeltaFHz', 428);

% -------------------- parse & prep --------------------
p = inputParser;
p.addParameter('DeltaFHz', 428);
p.addParameter('RefineB0', true);
p.addParameter('DoR2s', false);
p.addParameter('B0RangeHz', 300);
p.addParameter('NumIter', 3);
p.addParameter('MedianWin', 3);
p.parse(varargin{:});
deltaF = p.Results.DeltaFHz;
doRefine = p.Results.RefineB0;
doR2s = p.Results.DoR2s;
b0Range = p.Results.B0RangeHz;
nit = p.Results.NumIter;
medw = p.Results.MedianWin;

S = double(S);
sz = size(S);
if numel(TEs) < 3
    error('Need at least 3 echo times.');
end
TEs = TEs(:).';
S = S(:,:,1:3);
TEs = TEs(1:3);

[Ny,Nx,~] = size(S);
W = zeros(Ny,Nx);
F = zeros(Ny,Nx);
fB0 = zeros(Ny,Nx);
R2s = zeros(Ny,Nx);

% small helper
    function [Wc,Fc] = solveWF(Sv, TEs, f0, R2)
        % Given f0 and R2*, solve linear least-squares for complex W, F
        % Sv: [3x1] complex voxel data
        E0 = exp( 1i*2*pi*f0*TEs - R2*TEs );              % [1x3]
        A = [E0.' , (E0.*exp(1i*2*pi*deltaF*TEs)).'];     % [3x2]
        X = A \ Sv;                                       % 2x1
        Wc = X(1); Fc = X(2);
    end

% -------------------- initial fB0 estimate --------------------
% Alias-aware slope between echo phases; project into Nyquist of ΔTE
S1 = S(:,:,1); S2 = S(:,:,2); S3 = S(:,:,3);
dphi12 = angle(conj(S1).*S2);
dphi23 = angle(conj(S2).*S3);

dTE12 = TEs(2)-TEs(1);
dTE23 = TEs(3)-TEs(2);

% crude slope (average of two)
f_est = (dphi12/(2*pi)/dTE12 + dphi23/(2*pi)/dTE23)/2;

% clip to reasonable range to reduce wild outliers
f_est = max(min(f_est, b0Range), -b0Range);

% light spatial median filter to stabilize
if medw >= 3
    try
        f_est = medfilt2(f_est, [medw medw], 'symmetric');
    catch
        % if Signal Proc toolbox missing, skip filtering
    end
end

fB0 = f_est;

% -------------------- optional R2* one-shot (very mild) --------------------
if doR2s
    % Using log-magnitude linear fit after demodulating B0 (still approximate)
    Em = abs( S .* exp(-1i*2*pi*fB0.*reshape(TEs,1,1,[])) );
    y = log(max(Em, eps));
    X = [ones(Ny,Nx,1), reshape(TEs,1,1,[])];
    % least squares per-voxel: slope ~ -R2*
    % Solve small 3-point LS analytically:
    t1 = TEs(1); t2 = TEs(2); t3 = TEs(3);
    % Pre-compute normal equations constants
    A11 = 3;
    A12 = t1 + t2 + t3;
    A22 = t1^2 + t2^2 + t3^2;
    Det = A11*A22 - A12*A12;
    % sums
    sy  = sum(y,3);
    sty = t1*y(:,:,1) + t2*y(:,:,2) + t3*y(:,:,3);
    a0 = ( A22.*sy - A12.*sty ) ./ max(Det,eps);
    a1 = ( A11.*sty - A12.*sy ) ./ max(Det,eps);
    R2s = max(0, -a1);  % clamp negative to zero
else
    R2s(:) = 0;
end

% -------------------- alternating solve + B0 refinement --------------------
for it = 1:max(1, nit)
    % Solve W,F with current fB0,R2*
    for y = 1:Ny
        Sv = squeeze(S(y,:,:)).';
        for x = 1:Nx
            s = Sv(:,x); % [3x1]
            if all(s==0)
                continue;
            end
            [W(y,x),F(y,x)] = solveWF(s, TEs, fB0(y,x), R2s(y,x));
        end
    end

    if ~doRefine
        break;
    end

    % 1-D Gauss-Newton step on fB0: update that minimizes residual norm
    % r = S - (W + F*e^{i2πΔf TE}) * e^{i2π fB0 TE} * e^{-R2* TE}
    % Jacobian wrt fB0: J = -i 2π TE .* M , where M = model without error
    twoPi = 2*pi;
    for y = 1:Ny
        for x = 1:Nx
            if W(y,x)==0 && F(y,x)==0, continue; end
            E0 = exp( 1i*twoPi*fB0(y,x)*TEs - R2s(y,x)*TEs );
            M  = ( W(y,x) + F(y,x)*exp(1i*twoPi*deltaF*TEs) ).* E0; % [1x3]
            r  = squeeze(S(y,x,:)).' - M;                             % [1x3]
            J  = -1i*twoPi*TEs .* M;                                  % [1x3]
            % GN step: df = real((J^H r) / (J^H J))
            num = sum(conj(J).*r);
            den = sum(conj(J).*J) + 1e-12;
            df  = real(num/den);
            % damp small
            fB0(y,x) = fB0(y,x) + 0.7*df;
        end
    end

    if medw >= 3
        try
            fB0 = medfilt2(fB0, [medw medw], 'symmetric');
        catch
        end
    end
end

% -------------------- done --------------------
end
