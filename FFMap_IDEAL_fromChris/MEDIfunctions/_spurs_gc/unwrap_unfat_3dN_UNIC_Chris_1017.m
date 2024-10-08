%% Note: this is a new version for alternative jumpe move
% in the previous version we use differenct jump length for odd and even iter, which may bring invalid label

% in this new version: For nodes with fat jump by k+1,m-1
% For nodes without fat jump by m+1

% consider magnitude weighting 

%% Alternative jump move for joint phase unwrapping and unfat
% consider the case when ff = w*pi
% optimize y = x + 2k*pi + w*m*pi = x + (2k+wm) pi = x + K * pi
% optimize K, then compute k and m 
% use alternative jump move
% K(k,m)   0(0,0) w(0,1) 4(2,0) 4+w(2,1) 6(3,0) 6+w(3,1)


function [kappa,m_fat,unwph,iter,erglist,kappa_iter] = unwrap_unfat_3dN_UNIC_Chris_1017(voxel_size,iMag,w,psi,p,Mask,varargin)
%puma_ho   Graph optimization based phase unwrapping algorithm.
%   [unwph,iter,erglist] = puma_ho(psi,p,'potential',potential,'cliques',cliques, 'qualitymaps', qualitymaps,
%   'schedule',schedule)
%   Unwrapps the observed modulo-2pi phase, casting the problem as an energy minimization via graph mincut
%   calculation. The algorithm is described in "Phase Unwrapping via Graph Cuts" submitted to IEEE IP, October, 2005.
%   Herein we make a generalization, by allowing cliques of any order (not just 1st order).
%
%   Authors: Jose Bioucas-Dias and Gon�alo Valad�o
%
%   Last change: Goncalo Valadao (19/9/2012 23h46m)
%
% ======================== REQUIRED INPUT PARAMETERS ======================
% Parameter             Values
% name                  and description
% =========================================================================
%
% psi                   (double) The wrapped phase image.
% p                     (double) It defines the clique potential exponent (>0).
%
% ======================== OPTIONAL INPUT PARAMETERS ======================
% Parameter             Values
% name                  and description
% =========================================================================
%
% potential             (1x1 struct array) This struct array has 2 fields:
% potential.quantized   ('yes','no') Default: 'no'.
% potential.threshold   (double) it defines a region over which the
%                        potential grows quadratically. By default is pi.
%
% cliques               (nx2 double matrix) Each row defines the
%                       "displacement vector" corresponding to each clique.
%                       The first and the second columns correspond to
%                       cliques along rows and columns in the image, respecively.
%                       By default is [1 0;0 1] (first order cliques).
%
% qualitymaps           (size(psi) x n (nocliques) double array). The quality matrices
%                       may take values between 0 and 1 (value 1: discontinuity presence;
%                       value 0: discontinuity absence).
%                       There is one quality matrix per clique type. By default there is
%                       discontinuity absence. A quality value corresponding to a certain
%                       clique must be signalled in the pixel corresponding to the end of
%                       the clique displacement vector (for each pair of pixels).
%
% schedule              (double vector) This vector contains a schedule of jump sizes.
%
% verbose               ('yes', 'no') -> display the unwrapped  phase along  
%                        the iterations.
%                        Default = 'yes'
%
%
% ========================= OUTPUT PARAMETERS =============================
% Parameter             Values
% name                  and description
% =========================================================================
% unwph                 (double array) This is the unwrapped phase image.
% iter                  (double) This is the number of iterations the algorithm runs through
% erglist               (double vector) This is a sequence of the energies of the unwrapped
%                       phases along the algorithm run.
%
% =========================== EXAMPLES ====================================
%       Note: the optional arguments must be provided in pairs string+value;
%       those pairs may be entered in any order.
%
%       potential.quantized = 'no'; potential.threshold = 0.5;
%       [unwph,iter,erglist] = puma_ho(psi,2,'potential',potential)
%
%       potential.quantized = 'yes'; potential.threshold = 2;
%       cliques = [1 0; 0 1; -1 1];
%       [unwph,iter,erglist] = puma_ho(psi,1,'cliques',cliques,'potential',potential)
%
%       potential.quantized = 'no';
%       potential.threshold = 0.1; cliques = [1 1];
%       qualitymaps = ones(size(psi,1),size(psi,2))
%       [unwph,iter,erglist] = puma_ho(psi,p,'potential',potential,'cliques',cliques,'qualitymaps',qualitymaps)

% ========================== REFERENCES ===================================
%   For reference see:
%   J. Bioucas-Dias and G. Valad�o, "Phase Unwrapping via Graph Cuts"
%   IEEE Transactions Image Processing, 2007 (to appear).
%   The algorithm here coded corresponds to a generalization for any
%   cliques set (not only vertical and horizontal).
%
%   J. Bioucas-Dias and J. Leit�o, "The ZpiM Algorithm for Interferometric Image Reconstruction
%   in SAR/SAS", IEEE Transactions Image Processing, vol. 20, no. Y, 2001.
%
%   The technique here employed is also based on the article:
%   V. kolmogorov and R. Zabih, "What Energy Functions can be Minimized via Graph Cuts?",
%   European Conference on Computer Vision, May 2002.
% =========================================================================
%
%
%
% Modification: 
%     
%    1 - Fix a bug in the way discontinties were deal with 
%        (Gon�alo Valad�o, Sep.,2012)
%
%    2 - Change in the potential default parameters:
%        potential.quantized = 'no';
%        potential.threshold = pi;
%        (J Bioucas-Dias, Sep.,2012)
%
%    3 - Introdution of the verbose input parameter:
%        verbose = 'yes' -> display (iter, enerrg_actual, jump_size)
%                           and display the unwrapped  phase along
%                           the iterations.
%        Default = 'yes'
%        
%    
%   
%
% ------------------------------------------------------------------
% Author: Gon�alo Valad�o & Jose Bioucas-Dias, 2007
%
%

%
%% -------------------------------------------------------------------------
%
% Copyright (July, 2007):        Gon�alo Valad�o (gvaladao@lx.it.pt) 
%                                Jos?Bioucas-Dias (bioucas@lx.it.pt)
%
% PUMA_HO is distributed under the terms of
% the GNU General Public License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ---------------------------------------------------------------------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Default values
potential.quantized     = 'no';
potential.threshold     = pi;
cliques                 = [1 0; 0 1];
qualitymaps             = repmat(zeros(size(psi,1),size(psi,2)),[1,1,2,size(psi,3)]); 
qualitymaps_z = zeros(size(psi,1),size(psi,2),size(psi,3));

qual=0;
schedule                = 1;
verbose                 = 'yes';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test for number of required parameters

% Error out if there are not at least the two required input arguments

if (nargin-length(varargin)) ~= 6
    error('Wrong number of required parameters');
end

% Read the optional parameters
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
elseif length(varargin)~=0
    for i=1:2:(length(varargin)-1)
        % change the value of parameter
        switch varargin{i}
            case 'potential'        % potential definition
                potential = varargin{i+1};
            case 'cliques'          % cliques to consider
                cliques = varargin{i+1};
            case 'qualitymaps'      % the quality maps
                qualitymaps = varargin{i+1};
                qual = 1;
            case 'schedule'         % jump size schedule
                schedule = varargin{i+1};
           case 'verbose'           % display partial unwrappings
                verbose = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end
    end
end

if (qual==1)&&(size(qualitymaps,3)~=size(cliques,1))
    error('qualitymaps must be a 3D matrix whos 3D size is equal to no. cliques. Each plane on qualitymaps corresponds to a clique.');
end

% INPUT AND INITIALIZING %%%%%%%%%%%%%%%%%%%%%%%%
th    = getfield(potential,'threshold');
quant = getfield(potential,'quantized');


[m,n,h] = size(psi); % Size of input
kappa = zeros(m,n,h); % Initial labeling  when w is positive, initial labeling is zeros.
% 10.2019 C.H change of initial point
kappa_fat = mod(kappa,1)==1;
kappa_aux = kappa;
iter = 0;
erglist = [];

[cliquesm,cliquesn] = size(cliques); % Size of input cliques
if qual ==0
    qualitymaps             = repmat(zeros(size(psi,1),size(psi,2)),[1,1,size(cliques,1),size(psi,3)]);
end
disc_bar = 1 - qualitymaps;
disc_bar_z = 1 - qualitymaps_z; 


% "maxdesl" is the maximum clique length used.
maxdesl = max(max(abs(cliques)));
% We define "base" which is a mask having ones in the region of interest(psi) and zeros upon a passe-partout
% having a constant length maxdesl+1.

base_xy = zeros(2*maxdesl+2+m,2*maxdesl+2+n,h); 
for hh = 1:h
    base_xy(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1,hh) = ones(m,n);
end
base_z = zeros(m,n,h);
for hh = 1:h
    base_z(:,:,hh) = ones(m,n); 
end

      %% Temporary movie maker
%        mov = avifile('example5.avi');
        %% Temporary movie maker

% PROCESSING   %%%%%%%%%%%%%%%%%%%%%%%%
for jump_size = schedule
    possible_improvment = 1;
    % note: here kappa refer to K in our case,kappa can be non interger
%     erg_previous = energy_ho_3d_pi(kappa,psi,base_xy,base_z,p,cliques,disc_bar,disc_bar_z,th,quant);
    psi_Mask = psi;
    kappa_Mask = kappa;
    erg_previous = energy_ho_3dW2_pi_Chris(voxel_size,iMag,kappa_Mask,psi_Mask,base_xy,base_z,p,cliques,disc_bar,disc_bar_z,th,quant,Mask);
    
    %     erg_previous = 1e10;

    while possible_improvment
        
        iter = iter + 1;
        fprintf('1-jump move, iter %d \n', iter);
        
%         if mod(iter,2) == 1
%             binaryEdge_term = w*pi;
%         else
%             binaryEdge_term = (2-w)*pi;           
%         end
        
        erglist = [erglist erg_previous];
        remain = [];
        % Here we put a passe-partout (constant length = maxdesl+1) in the images kappa and psi
%         base_kappa = zeros(2*maxdesl+2+m,2*maxdesl+2+n); 
%         base_kappa(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1) = kappa;
        
        base_kappa_xy = zeros(2*maxdesl+2+m,2*maxdesl+2+n,h);
        % the 0,1 matrix to note wheter one voxel has fat
        base_kappa_fat_xy = zeros(2*maxdesl+2+m,2*maxdesl+2+n,h); 
        
        for hh = 1:h
            base_kappa_xy(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1,hh) = kappa(:,:,hh);
            base_kappa_fat_xy(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1,hh) = kappa_fat(:,:,hh);
        end
        base_kappa_jump_xy = (2-w)*base_kappa_fat_xy + w*(1-base_kappa_fat_xy);
        
        base_kappa_z = kappa;
        
        
                
%         psi_base = zeros(2*maxdesl+2+m,2*maxdesl+2+n); 
%         psi_base(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1) = psi;
        psi_base_xy      = zeros(2*maxdesl+2+m,2*maxdesl+2+n,h); 
        for hh = 1:h
            psi_base_xy(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1,hh) = psi(:,:,hh);
        end 
        
        iMag_base = zeros(2*maxdesl+2+m,2*maxdesl+2+n,h); 
        for hh = 1:h
            iMag_base(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1,hh) = iMag(:,:,hh);
        end 

        
        %%%%% Added by Goncalo Valadao on 19/11/2012 22:15 %%%%%
        z = size(disc_bar,3);
        base_disc_bar  = repmat(zeros(2*maxdesl+2+m,2*maxdesl+2+n),[1 1 z size(psi,3)]); 
        base_disc_bar(maxdesl+2:maxdesl+2+m-1,maxdesl+2:maxdesl+2+n-1,:,:) = disc_bar;
        %%%%% Added by Goncalo Valadao on 19/11/2012 22:15 %%%%%
        
        for zz = 1:size(psi,3)
             for t = 1:cliquesm
                % The allowed start and end pixels of the "interpixel" directed edge
                base_start(:,:,t,zz) = circshift(base_xy(:,:,zz),[-cliques(t,1),-cliques(t,2)]).*base_xy(:,:,zz);
                base_end(:,:,t,zz) = circshift(base_xy(:,:,zz),[cliques(t,1),cliques(t,2)]).*base_xy(:,:,zz);

                % By convention the difference images have the same size as the
                % original ones; the difference information is retrieved in the
                % pixel of the image that is subtracted (end of the diff vector)
                auxili(:,:,zz) = circshift(base_kappa_xy(:,:,zz),[cliques(t,1),cliques(t,2)]);
                t_dkappa(:,:,t,zz) = (base_kappa_xy(:,:,zz)-auxili(:,:,zz));
                
                % the jump length of every voxels moved to smaller rows 
                auxili_jump(:,:,zz) = circshift(base_kappa_jump_xy(:,:,zz),[cliques(t,1),cliques(t,2)]);
                
                t_dkappa00(:,:,t,zz) = zeros(size(base_kappa_jump_xy(:,:,zz)));
                t_dkappa01(:,:,t,zz) = base_kappa_jump_xy(:,:,zz) - zeros(size(base_kappa_jump_xy(:,:,zz)));
                t_dkappa10(:,:,t,zz) = zeros(size(auxili_jump(:,:,zz))) - auxili_jump(:,:,zz);
                t_dkappa11(:,:,t,zz) = base_kappa_jump_xy(:,:,zz) - auxili_jump(:,:,zz);
                
                
                auxili2(:,:,zz) = circshift(psi_base_xy(:,:,zz),[cliques(t,1),cliques(t,2)]);
                dpsi(:,:,zz) = auxili2(:,:,zz) - psi_base_xy(:,:,zz);
                
                % add the weight of iMag
                auxiliW(:,:,zz) = circshift(iMag_base(:,:,zz),[cliques(t,1),cliques(t,2)]);
                diMag(:,:,t,zz) = (auxiliW(:,:,zz) + iMag_base(:,:,zz))/2;

                % Beyond base, we must multiply by
                % circshift(base,[cliques(t,1),cliques(t,2)]) in order to
                % account for frontier pixels that can't have links outside ROI

                %%%%%%%%  Changed on 19/11/2012 18:30 by Goncalo Valadao  %%%%%%%%

    %             a(:,:,t) = (2*pi*t_dkappa(:,:,t)-dpsi).*base.*circshift(base,[cliques(t,1),cliques(t,2)]);
    %             A(:,:,t) = clique_energy_ho(abs(a(:,:,t)),p,th,quant).*base.*circshift(base,[cliques(t,1),cliques(t,2)]);
    %             D(:,:,t) = A(:,:,t);
    %             C(:,:,t) = clique_energy_ho(abs(2*pi*jump_size + a(:,:,t)),p,th,quant).*base.*circshift(base,[cliques(t,1),cliques(t,2)]);
    %             B(:,:,t) = clique_energy_ho(abs(-2*pi*jump_size + a(:,:,t)),p,th,quant).*base.*circshift(base,[cliques(t,1),cliques(t,2)]);

%                 a(:,:,t,zz) = (1*pi*t_dkappa(:,:,t,zz)-dpsi(:,:,zz)).*base_xy(:,:,zz).*circshift(base_xy(:,:,zz),[cliques(t,1),cliques(t,2)]);
%                 A(:,:,t,zz) = clique_energy_ho(abs(a(:,:,t,zz)),p,th,quant).*base_xy(:,:,zz).*circshift(base_xy(:,:,zz),[cliques(t,1),cliques(t,2)]).*base_disc_bar(:,:,t,zz);
%                 D(:,:,t,zz) = A(:,:,t,zz);
%                 C(:,:,t,zz) = clique_energy_ho(abs(binaryEdge_term*jump_size + a(:,:,t,zz)),p,th,quant).*base_xy(:,:,zz).*circshift(base_xy(:,:,zz),[cliques(t,1),cliques(t,2)]).*base_disc_bar(:,:,t,zz);
%                 B(:,:,t,zz) = clique_energy_ho(abs(-binaryEdge_term*jump_size + a(:,:,t,zz)),p,th,quant).*base_xy(:,:,zz).*circshift(base_xy(:,:,zz),[cliques(t,1),cliques(t,2)]).*base_disc_bar(:,:,t,zz);
                
                a(:,:,t,zz) = (1*pi*t_dkappa(:,:,t,zz)-dpsi(:,:,zz)).*base_xy(:,:,zz).*circshift(base_xy(:,:,zz),[cliques(t,1),cliques(t,2)]); 
                b(:,:,t,zz) = (1*pi*(t_dkappa(:,:,t,zz) + t_dkappa01(:,:,t,zz))-dpsi(:,:,zz)).*base_xy(:,:,zz).*circshift(base_xy(:,:,zz),[cliques(t,1),cliques(t,2)]);
                c(:,:,t,zz) = (1*pi*(t_dkappa(:,:,t,zz) + t_dkappa10(:,:,t,zz))-dpsi(:,:,zz)).*base_xy(:,:,zz).*circshift(base_xy(:,:,zz),[cliques(t,1),cliques(t,2)]);
                d(:,:,t,zz) = (1*pi*(t_dkappa(:,:,t,zz) + t_dkappa11(:,:,t,zz))-dpsi(:,:,zz)).*base_xy(:,:,zz).*circshift(base_xy(:,:,zz),[cliques(t,1),cliques(t,2)]);
                
                A(:,:,t,zz) = diMag(:,:,t,zz).*clique_energy_ho(abs(a(:,:,t,zz)/voxel_size(1)),p,th,quant).*base_xy(:,:,zz).*circshift(base_xy(:,:,zz),[cliques(t,1),cliques(t,2)]).*base_disc_bar(:,:,t,zz);
                B(:,:,t,zz) = diMag(:,:,t,zz).*clique_energy_ho(abs(b(:,:,t,zz)/voxel_size(1)),p,th,quant).*base_xy(:,:,zz).*circshift(base_xy(:,:,zz),[cliques(t,1),cliques(t,2)]).*base_disc_bar(:,:,t,zz);
                C(:,:,t,zz) = diMag(:,:,t,zz).*clique_energy_ho(abs(c(:,:,t,zz)/voxel_size(1)),p,th,quant).*base_xy(:,:,zz).*circshift(base_xy(:,:,zz),[cliques(t,1),cliques(t,2)]).*base_disc_bar(:,:,t,zz);
                D(:,:,t,zz) = diMag(:,:,t,zz).*clique_energy_ho(abs(d(:,:,t,zz)/voxel_size(1)),p,th,quant).*base_xy(:,:,zz).*circshift(base_xy(:,:,zz),[cliques(t,1),cliques(t,2)]).*base_disc_bar(:,:,t,zz);
                %%%%%%%%  End of changes on 19/11/2012 18:30 by Goncalo Valadao  %%%%%%%%  

                % The circshift by [-cliques(t,1),-cliques(t,2)] is due to the fact that differences are retrieved in the
                % "second=end" pixel. Both "start" and "end" pixels can have source and sink connections.
                source(:,:,t,zz) = circshift((C(:,:,t,zz)-A(:,:,t,zz)).*((C(:,:,t,zz)-A(:,:,t,zz))>0),[-cliques(t,1),-cliques(t,2)]).*base_start(:,:,t,zz);
                sink(:,:,t,zz) = circshift((A(:,:,t,zz)-C(:,:,t,zz)).*((A(:,:,t,zz)-C(:,:,t,zz))>0),[-cliques(t,1),-cliques(t,2)]).*base_start(:,:,t,zz);

                source(:,:,t,zz) = source(:,:,t,zz) + ((D(:,:,t,zz)-C(:,:,t,zz)).*((D(:,:,t,zz)-C(:,:,t,zz))>0)).*base_end(:,:,t,zz);
                sink(:,:,t,zz) = sink(:,:,t,zz) + ((C(:,:,t,zz)-D(:,:,t,zz)).*((C(:,:,t,zz)-D(:,:,t,zz))>0)).*base_end(:,:,t,zz);
            end
        end
        
        

        % We get rid of the "pass-partous"
%         for zz = 1:size(psi,3)
%             source(1:maxdesl+1,:,:,zz)=[]; source(m+1:m+maxdesl+1,:,:,zz)=[]; source(:,1:maxdesl+1,:,zz)=[]; source(:,n+1:n+maxdesl+1,:,zz)=[];
%             sink(1:maxdesl+1,:,:,zz)=[]; sink(m+1:m+maxdesl+1,:,:,zz)=[];sink(:,1:maxdesl+1,:,zz)=[]; sink(:,n+1:n+maxdesl+1,:,zz)=[];
%         end
        
        for zz = 1:size(psi,3)
            source2(:,:,:,zz) = source(maxdesl+2:m+maxdesl+1,maxdesl+2:n+maxdesl+1,:,zz);
            sink2(:,:,:,zz) = sink(maxdesl+2:m+maxdesl+1,maxdesl+2:n+maxdesl+1,:,zz);
        end
        source = source2;
        sink = sink2;
                        
        auxiliar1 = B + C - A - D;
        non_sub = sum(auxiliar1(:)<0);
        
        for zz = 1:size(psi,3)
%             auxiliar1(1:maxdesl+1,:,:,zz)=[]; auxiliar1(m+1:m+maxdesl+1,:,:,zz)=[]; auxiliar1(:,1:maxdesl+1,:,zz)=[]; auxiliar1(:,n+1:n+maxdesl+1,:,zz)=[];
%             base_start(1:maxdesl+1,:,:,zz)=[]; base_start(m+1:m+maxdesl+1,:,:,zz)=[]; base_start(:,1:maxdesl+1,:,zz)=[]; base_start(:,n+1:n+maxdesl+1,:,zz)=[];
%             base_end(1:maxdesl+1,:,:,zz)=[]; base_end(m+1:m+maxdesl+1,:,:,zz)=[]; base_end(:,1:maxdesl+1,:,zz)=[]; base_end(:,n+1:n+maxdesl+1,:,zz)=[];
            auxiliar1_2(:,:,:,zz) = auxiliar1(maxdesl+2:m+maxdesl+1,maxdesl+2:n+maxdesl+1,:,zz);
            base_start2(:,:,:,zz) = base_start(maxdesl+2:m+maxdesl+1,maxdesl+2:n+maxdesl+1,:,zz);
            base_end2(:,:,:,zz) = base_end(maxdesl+2:m+maxdesl+1,maxdesl+2:n+maxdesl+1,:,zz);
        end
        auxiliar1 = auxiliar1_2;
        base_start = base_start2;
        base_end = base_end2;
        % We construct the "remain" and the "sourcesink" matrices
        
        for zz = 1:size(psi,3)
            for t=1:cliquesm
                start = find(base_start(:,:,t,zz)~=0); 
                endd = find(base_end(:,:,t,zz)~=0);
                auxiliar2 = auxiliar1(:,:,t,zz);
                auxiliar3 = [start+(zz-1)*m*n endd+(zz-1)*m*n  auxiliar2(endd).*(auxiliar2(endd)>0) zeros(size(endd,1),1)];
                %auxiliar3 = [start endd  auxiliar2(endd) zeros(size(endd,1),1)];
                remain = [remain; auxiliar3];
            end
            
        end 
        
        
        %remain = sortrows(remain,[1 2]);
%         sourcefinal = sum(source,3);
%         sinkfinal = sum(sink,3);
%         sourcesink = [(1:m*n)' sourcefinal(:) sinkfinal(:)];
        sourcesink = [];
        for zz = 1:size(psi,3)
            sourcefinal(:,:,zz) = sum(source(:,:,:,zz),3);
            sinkfinal(:,:,zz) = sum(sink(:,:,:,zz),3);
            s_temp = sourcefinal(:,:,zz);
            t_temp = sinkfinal(:,:,zz);
            sourcesink = [sourcesink; (1:m*n)'+m*n*(zz-1) s_temp(:) t_temp(:)];
        end
        
        % TO-DO: add the energy terms of edges in the z direction
        base_disc_bar_z = disc_bar_z;
        
        base_kappa_z = kappa;
        % fat label and the jump difference 
        base_kappa_fat_z = kappa_fat;
      
        base_kappa_jump_z = (2-w)*base_kappa_fat_z + w*(1-base_kappa_fat_z);
        
        % find one mistake! the following should be in the for loop
%         t_dkappa00z(:,:,zz) = zeros(size(base_kappa_jump_z(:,:,zz)));
%         t_dkappa01z(:,:,zz) = base_kappa_jump_z(:,:,zz) - zeros(size(base_kappa_jump_z(:,:,zz-1)));
%         t_dkappa10z(:,:,zz) = zeros(size(base_kappa_jump_z(:,:,zz))) - base_kappa_jump_z(:,:,zz-1);
%         t_dkappa11z(:,:,zz) = base_kappa_jump_z(:,:,zz) - base_kappa_jump_z(:,:,zz-1);
        
        
        dkappa(:,:,1) = zeros(m,n);
        dkappa_fat(:,:,1) = zeros(m,n);
        dpsi_z(:,:,1) = zeros(m,n);
        
         % add the weight of iMag
        diMagz(:,:,1) = zeros(m,n);
        
        for zz = 2:size(psi,3) 
%             diMagz(:,:,zz) = ((1./iMag(:,:,zz-1)).^2 + (1./iMag(:,:,zz)).^2).^(-0.5);
            diMagz(:,:,zz) = (iMag(:,:,zz-1) + iMag(:,:,zz))/2;
        end
        
        AA(:,:,1) = zeros(m,n);
        BB(:,:,1) = zeros(m,n);
        CC(:,:,1) = zeros(m,n);
        DD(:,:,1) = zeros(m,n);
        source_z = zeros(m,n,h);
        sink_z = zeros(m,n,h);
        for zz = 2:size(psi,3)
            % the plane with higher z minus the plane with one smaller plane
            dkappa(:,:,zz) = base_kappa_z(:,:,zz) - base_kappa_z(:,:,zz-1);
            % the voxel potential in smaller z minus the ones with one
            % bigger ones
            dpsi_z(:,:,zz) = psi(:,:,zz-1) - psi(:,:,zz);   
            
            t_dkappa00z(:,:,zz) = zeros(size(base_kappa_jump_z(:,:,zz)));
            t_dkappa01z(:,:,zz) = base_kappa_jump_z(:,:,zz) - zeros(size(base_kappa_jump_z(:,:,zz-1)));
            t_dkappa10z(:,:,zz) = zeros(size(base_kappa_jump_z(:,:,zz))) - base_kappa_jump_z(:,:,zz-1);
            t_dkappa11z(:,:,zz) = base_kappa_jump_z(:,:,zz) - base_kappa_jump_z(:,:,zz-1);
           
            aa(:,:,zz) = (1*pi*(dkappa(:,:,zz) + t_dkappa00z(:,:,zz))-dpsi_z(:,:,zz)).*base_disc_bar_z(:,:,zz);
            bb(:,:,zz) = (1*pi*(dkappa(:,:,zz) + t_dkappa01z(:,:,zz))-dpsi_z(:,:,zz)).*base_disc_bar_z(:,:,zz);
            cc(:,:,zz) = (1*pi*(dkappa(:,:,zz) + t_dkappa10z(:,:,zz))-dpsi_z(:,:,zz)).*base_disc_bar_z(:,:,zz);
            dd(:,:,zz) = (1*pi*(dkappa(:,:,zz) + t_dkappa11z(:,:,zz))-dpsi_z(:,:,zz)).*base_disc_bar_z(:,:,zz);
            % error found:  the * should be .*        AA(:,:,zz) =
            % clique_energy_ho(abs(b(:,:,zz)),p,th,quant).*base_z(:,:,zz)*base_disc_bar_z(:,:,zz);
            AA(:,:,zz) = diMagz(:,:,zz).*clique_energy_ho(abs(aa(:,:,zz)/voxel_size(3)),p,th,quant).*base_z(:,:,zz).*base_disc_bar_z(:,:,zz);
            BB(:,:,zz) = diMagz(:,:,zz).*clique_energy_ho(abs(bb(:,:,zz)/voxel_size(3)),p,th,quant).*base_z(:,:,zz).*base_disc_bar_z(:,:,zz);
            CC(:,:,zz) = diMagz(:,:,zz).*clique_energy_ho(abs(cc(:,:,zz)/voxel_size(3)),p,th,quant).*base_z(:,:,zz).*base_disc_bar_z(:,:,zz);
            DD(:,:,zz) = diMagz(:,:,zz).*clique_energy_ho(abs(dd(:,:,zz)/voxel_size(3)),p,th,quant).*base_z(:,:,zz).*base_disc_bar_z(:,:,zz);
%             DD(:,:,zz) = AA(:,:,zz);
%             CC(:,:,zz) = clique_energy_ho(abs(binaryEdge_term*jump_size + b(:,:,zz)),p,th,quant).*base_z(:,:,zz).*base_disc_bar_z(:,:,zz);
%             BB(:,:,zz) = clique_energy_ho(abs(-binaryEdge_term*jump_size + b(:,:,zz)),p,th,quant).*base_z(:,:,zz).*base_disc_bar_z(:,:,zz);
            
            source_z(:,:,zz-1) = source_z(:,:,zz-1) + (CC(:,:,zz)-AA(:,:,zz)).*((CC(:,:,zz)-AA(:,:,zz))>0);
            sink_z(:,:,zz-1) = sink_z(:,:,zz-1) + (AA(:,:,zz)-CC(:,:,zz)).*((AA(:,:,zz)-CC(:,:,zz))>0);

            source_z(:,:,zz) = source_z(:,:,zz) + ((DD(:,:,zz)-CC(:,:,zz)).*((DD(:,:,zz)-CC(:,:,zz))>0));
            sink_z(:,:,zz) = sink_z(:,:,zz) + ((CC(:,:,zz)-DD(:,:,zz)).*((CC(:,:,zz)-DD(:,:,zz))>0));
        end
        sourcesink2 = [];
        for zz = 1:size(psi,3)
           
            s_temp = source_z(:,:,zz);
            t_temp = sink_z(:,:,zz);
            sourcesink2 = [sourcesink2; (1:m*n)'+m*n*(zz-1) s_temp(:) t_temp(:)];
        end
        
        sourcesink = sourcesink + sourcesink2;
        sourcesink(:,1) =  sourcesink(:,1)/2; 
        
         
        
        
        auxiliar1 = BB + CC - AA - DD;
        non_sub = sum(auxiliar1(:)<0);
        
        base = ones(m,n);     
        for zz = 2:size(psi,3)
            start = find(base~=0) + (zz-2)*m*n;
            endd = find(base~=0) + (zz-1)*m*n;
            auxiliar2 = auxiliar1(:,:,zz);         
            auxiliar3 = [start endd  auxiliar2(:).*(auxiliar2(:)>0) zeros(size(endd,1),1)];
            %auxiliar3 = [start endd  auxiliar2(endd) zeros(size(endd,1),1)];
            remain = [remain; auxiliar3];
        end
        

        % KAPPA RELABELING
        [flow,cutside] = mincut(sourcesink,remain);
        % The relabeling of each pixel relies on the cutside of that pixel:
        %       if the cutside = 0 (source) then we increment the label by jump_size
        %       if the cutside = 1 (sink) then the label remains unchanged
        Majorpix=sum(cutside(Mask,2))/sum(Mask(:));%major pixel
        kappa_temp = zeros(size(kappa));
        kappa_temp(cutside(:,1)) = cutside(:,2);
        %if Majorpix>0.5 % major pixel are 1-> the zero part should jump +
        %    kappa_temp(cutside(:,1)) = kappa(cutside(:,1)) + (1 - cutside(:,2))*jump_size;
        %else %major pixels are 0 -> the one part should jump -
        %    kappa_temp(cutside(:,1)) = kappa(cutside(:,1)) - cutside(:,2)*jump_size;
        %end
        % C.H 2019/10/17 sign corrction
        %kappa_temp = 1-kappa_temp;
        kappa_aux = kappa + kappa_temp.*((2-w)*kappa_fat + w*(1-kappa_fat));
        
%          kappa_temp(kappa_binary == 1) = 1;
%         if mod(iter,2) ==1
%             % cutside is 0, move w; otherwise, do not move
%             kappa_aux(cutside(:,1)) = kappa(cutside(:,1)) + (1 - cutside(:,2))*jump_size*w;
%             
%         else
% %           kappa_aux(cutside(:,1)) = kappa(cutside(:,1)) + 2*(1 - cutside(:,2))*jump_size + cutside(:,2)*jump_size*w;
%             % cutside is 0, move (2-w); otherwise, do not move
%             kappa_aux(cutside(:,1)) = kappa(cutside(:,1)) + 2*(1 - cutside(:,2))*jump_size - (1 - cutside(:,2))*jump_size*w;
%             
%         end
%         kappa_aux(cutside(:,1)) = kappa(cutside(:,1)) + (1 - cutside(:,2))*jump_size;
        % CHECK ENERGY IMPROVEMENT
%       erg_actual = energy_ho(kappa_aux,psi,base,p,cliques,disc_bar,th,quant);
%       erg_actual = energy_ho_3d_pi(kappa_aux,psi,base_xy,base_z,p,cliques,disc_bar,disc_bar_z,th,quant);
        psi_Mask = psi;
        kappa_aux_Mask = kappa_aux; 
        erg_actual = energy_ho_3dW2_pi_Chris(voxel_size,iMag,kappa_aux_Mask,psi_Mask,base_xy,base_z,p,cliques,disc_bar,disc_bar_z,th,quant,Mask);
        %
        if (erg_actual < erg_previous)
            erg_previous = erg_actual;
            kappa = kappa_aux.*Mask;
            kappa_iter(:,:,:,iter) = kappa;
            kappa_fat = mod(kappa,1)== 1;
        else
           
            possible_improvment = 0;
            kappa_iter(:,:,:,iter) = kappa;
            unwph = 1*pi*kappa.*Mask + psi;
          
            m_fat = 1 - (mod(kappa,1)==0);
            
%             kappa(mod(kappa,1)==1) ;
%             m_fat = mod(kappa,2);
%             k_wrap = (kappa - m_fat)/2;

        end
%         mesh( 2*pi*kappa + psi)
%         view(-30,30);kappa
%         surfl(2*pi*kappa + psi);shading interp; colormap(gray);
        
%         if strcmp(verbose,'yes')
%            imagesc(2*pi*kappa + psi); colormap(gray);
%            % display to the current figure
%            figure(gcf);
%            drawnow;
%         end
        

        %surfl(2*pi*kappa + psi);shading interp; 
        %colormap(gray);
        %imagesc( 2*pi*kappa + psi);
        

        
        
        %% Temporary movie maker
%        F=getframe(gca);
%        mov=addframe(mov,F);
        %% Temporary movie maker        
        
        %drawnow;
        %iter %#ok<NOPRT>
        %erg_actual %#ok<NOPRT>
        %jump_size %#ok<NOPRT>
        clear base_start base_end source sink auxiliar1 auxiliar2 A B C D AA BB CC DD;
    end % while
%     title('Puma solution');
end %for

% Temporary movie maker
%mov=close(mov);
% Temporary movie maker
