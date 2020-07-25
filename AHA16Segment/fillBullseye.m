% SURFACEOBJECT = FILLBULLSEYE(varargin)  creates a 
% wedge "patch" using xdata/ydata from a polar plot. The color of the patch 
% is determined scaled by the value of the vector. To be used with
% createBullseye.
% 
% 
% varargin can take the following inputs:
% -cdata,rho1,rho2,theta1,theta2
% -cdata (assumes 
% 
% ROWS of cdata should correspond to the same "rho" in the bullseye. The
%   first row represents the outermost bullseye ring.
% COLS of cdata should correspond to the "theta" in the bullseye.
% 
% =========================================================================
% Adrian Lam                                                  Oshinski Lab
% August 7 2014                   
% =========================================================================


function surfaceObject = fillBullseye(varargin)

    cdata = varargin{1};
    
    if iscolumn(cdata)
        %Vector needs to be in columns for upsampling to work properly.
        cdata = cdata';
    end
            
    sz = size(cdata);
    
    if nargin == 1
        
        % Defaults are made to fit the AHA bullseye. 
        rho1 = 0.5;
        rho2 = 2;
        theta1 = 0;
        theta2 = 360;        

   elseif nargin == 5
        
        rho1 = varargin{2};
        rho2 = varargin{3};
        theta1 = varargin{4};
        theta2 = varargin{5};
        
    end
    
    dtheta = (theta2 - theta1)/sz(2);
    
    % dtheta should be less than 5 degrees for a smooth circle
    upsmp = 1;
    while dtheta > 5
        
        upsmp = upsmp + 1;
        dtheta = (theta2 - theta1) / ( sz(2) * upsmp );        
        
    end
    
    %Upsample cdata to match dtheta. No interpolation performed, but can be
    %added.
    cdata = reshape(repmat(cdata,upsmp,1),sz(1),sz(2)*upsmp);    
        
    drho = (rho2 - rho1)/sz(1);
    theta = repmat(theta1:dtheta:theta2,sz(1)+1,1);
    lin = rho2:-drho:rho1;
    
    % Create meshgrid for surface plot.
    X = repmat(lin',1,sz(2)*upsmp+1) .* cos(theta*pi/180);
    Y = repmat(lin',1,sz(2)*upsmp+1) .* sin(theta*pi/180);
        
    surfaceObject = surf(gca,X,Y,zeros(size(X)),cdata);
    set(surfaceObject,'EdgeColor','none');
    set(gca,'View',[0 90]);

end