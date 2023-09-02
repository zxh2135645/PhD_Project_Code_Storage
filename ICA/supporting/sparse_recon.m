U = U(:);

morozov=msdev^2*numel(kspace_data)

mov = @(x)reshape(x,Ny,Nx,[]);
% W=@(x)dispim(mov(reshape(x,[],L)*Wt));
% Wh=@(x)vec(reshape(padarray(x,[ovs/2, ovs/2, 0]),[],size(Wt,2))*Wt');
% WhW=@(x)vec(padarray(dispim(mov(reshape(x,[],L)*WtWh)),[ovs/2, ovs/2, 0]));
W=@(x)mov(reshape(x,[],L)*Wt);
Wh=@(x)vec(reshape(x,[],size(Wt,2))*Wt');
WhW=@(x)vec(mov(reshape(x,[],L)*WtWh));
figmake=@(WU)sum(abs(WU),3);

WU=W(U);
lambda=2*msdev^2/mean(abs(WU(:)-median(real(WU(:)))-1i*median(imag(WU(:)))))
Y=zeros(size(WU));
Z = WU;
alpha=max(abs(WU(:)))
figure,hist(abs(WU(:)),1000);
im0 = log(figmake(WU));
figure,imshow(im0,[]);
fprintf('Adjust lambda and alpha if desired, then return.\n')
keyboard;
rho=lambda/alpha;
tic;
for it=2:25
    WU = W(U);
    figure(100),imshow([im0, log(figmake(WU))],[]),drawnow;
    figure(101),imshow(5*[exp(im0),figmake(WU)]/max(exp(im0(:)))),drawnow;
    %   Zold = Z;
    Z=WU+Y/rho;
    Z = sign(Z).*max(abs(Z)-alpha,0);
    Y=Y+rho*(WU-Z);
    %   pr=norm(WU(:)-Z(:));
    %   dr=rho*norm(vec(Wh(Z-Zold)));
    %   step=pr/dr
    clear Zold;
    
    %   if step < 0.01
    %       step = 1
    %   else
    %       step = median([1.25 step 2])
    %   end
    step = 1.25;
    alpha=alpha/step;
    rho=lambda/alpha;
    
    Uold=U;
    switch Trajectory
        case 'Radial'
            U=pcg(@(x)AhA_ps(x,st,SEs,Phi_rt,sp)+rho/2*WhW(x),...
                Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 20]),M,[],U);
        case 'Linogram'
            U=pcg(@(x)AhA_ps_lin(x,st,SEs,Phi_rt,A,Ah,thetas)+rho/2*WhW(x),...
                Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 20]),M,[],U);
        case 'Cartesian'
            U=pcg(@(x)AhA_ps_cart(x,st,SEs,Phi_rt,Ny,Nx,Nfft,shifter)+rho/2*WhW(x),...
                Ahb(:)+rho/2*Wh(Z-Y/rho),[],median([5 it 20]),M,[],U);
    end
    eps = norm(U(:)-Uold(:))/norm(Uold(:))
    clear Uold;
    if (eps < 1e-3) && (eps ~= 0)
        break;
    end
end
toc;