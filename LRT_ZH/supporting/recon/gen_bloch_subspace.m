Nseg = params.lSegments;
% if ~strcmp(ScanType,'Cine')
%   alpha_deg = params.adFlipAngleDegree;
%   switch ScanType
%     case 'IR'
%       alpha0_deg = 180;
%       cL = 5;
%     case 'T2IR'
%       alpha0_deg='T2IR';
%       cL=5;
%     case 'T2prep'
%       alpha0_deg = 'T2';
%       cL = 5;
%     case {'SR',}
%       alpha0_deg = 90;
%       %%%%%%%%%%%%%%%%%%
%       cL = 3; % disp('forced cL=2, but cL=3 would be better (with 4 navs/SR block)')
%       %%%%%%%%%%%%%%%%%%
%   end
%   
%   [curveU,curveS]=gen_curve_subspace(Nseg,params.lEchoSpacing,alpha_deg,alpha0_deg,(params.alTR_seconds-params.lEchoSpacing*Nseg)/2);
%   
%   curvePhi = curveU(:,1:cL);
%   
% else
  cL = 1;
%   curvePhi = [];
% end