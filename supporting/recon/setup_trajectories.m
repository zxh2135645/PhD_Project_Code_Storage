switch Trajectory
  case 'Radial'
    fib = [233 377]; %hard coded in CV_stars
    trajs = fib(2);
    
    % phi=(1+sqrt(5))/2; %exact
    phi = fib(2)/fib(1);
    theta=180/phi;
    
    %Find partition order
    ParOrder = twix_obj.image.Par(:);
    ParOrder = ParOrder((cutoff+1):cutoff_end,:,:,:);
    ParOrder(nav_indices)=[];
    
    thetas=mod(theta*mod(sum(mod(1:cutoff,SGblock)~=1)+(1:size(kspace_data,1)),2*trajs),360);
    
    %gradient delay correction
    cand_coords=unique([thetas.',ParOrder],'rows'); %set of (theta,kz)-coordinates
    cand_coords=intersect(cand_coords(cand_coords(:,1)<180,:),...
      bsxfun(@minus,cand_coords(cand_coords(:,1)>=180,:),[180 0]),'rows'); %set of coordinates which also have opposed pairings
    ksig=zeros(size(cand_coords,1),size(kspace_data,2),size(kspace_data,4));
    nsig=ksig;
    for j=1:size(cand_coords,1)
      t_ind=(thetas==(cand_coords(j,1)+180)) & (ParOrder.'==cand_coords(j,2));
      ksig(j,:,:)=mean(ifft(kspace_data(t_ind,end:-1:1,:),[],2),1); %rather than mean, rank-1 correction would be better
      t_ind=(thetas==cand_coords(j,1)) & (ParOrder.'==cand_coords(j,2));
      nsig(j,:,:)=mean(ifft(kspace_data(t_ind,:,:),[],2),1);
    end
    x=[0:(Norig-1) (-Norig):-1];
    cost=@(shift)norm(reshape(ksig-nsig.*repmat(exp(1i*(2*shift-pi/Norig)*x),[size(ksig,1) 1 Ncoils]),[],1));
    shifts=linspace(-.1,.1,201);
    mincost=inf;
    for j=1:numel(shifts)
      tempcost=cost(shifts(j));
      if tempcost < mincost
        mincost=tempcost;
        shift=shifts(j);
      end
    end
    shift=fminsearch(cost,shift);
    
    kspace_data=fft(ifft(kspace_data,[],2).*repmat(exp(1i*shift*x),...
      [size(kspace_data,1) 1 size(kspace_data,3) Ncoils]),[],2);
    nav_data=fft(ifft(nav_data,[],2).*repmat(exp(1i*shift*x),...
      [size(nav_data,1) 1 size(nav_data,3) Ncoils]),[],2);
    
    kspace_data(thetas>=180,:,:,:)=kspace_data(thetas>=180,[1 end:-1:2],:,:);
    thetas=mod(thetas,180);
    klast = floor(size(kspace_data,1)/trajs)*trajs;
    thetas(1:klast)=repmat(thetas(1:trajs),[1 klast/trajs]);
    thetas(klast+1:end)=thetas(1:end-klast);
    
    r=linspace(-pi,pi,2*Norig+1); r(end)=[];
    om1=cos(thetas(1:trajs)*pi/180).'*r;
    om2=sin(thetas(1:trajs)*pi/180).'*r;
    om=[om1(:), om2(:)];
    
    fprintf('Setting up NUFFT...')
    if useGPU
      st.M=size(om,1);
      st.om=om;
      st.F=[];
      st.Nd=[N N Nz];
      GhG=[];
      
    else
        st=nufft_init(om,[N N],[6 6],2*[N N],[N N]/2);
        GhG = st.p.G'*st.p.G;
    end
    clear om1 om2 om
    fprintf('done.\n');
    
    case 'Cartesian'
        pe_indices(:,1) = pe_indices(:,1) * 2;
        ParOrder = pe_indices(:,2);
        
        [tempy, tempz]=ndgrid(1:Ny,1:Nz);
        tempy=fftshift(tempy);
        tempz=fftshift(tempz);
        temp=sub2ind([Ny Nz],tempy,tempz);
        st.pes=sub2ind([Ny Nz],pe_indices(:,1), ParOrder(:));
        st.pes=temp(st.pes);
        
        st.w=zeros(Ny*Nz,1);
        st.Ahmat=spalloc(Ny*Nz,numel(st.pes),numel(st.pes));
        for j=1:Ny*Nz
            st.w(j)=sum(st.pes==j);
            st.Ahmat(j,:)=(st.pes.'==j);
        end
        
        st.w(1)=st.w(1)-size(nav_data,1);
        st.winv = 1./st.w;
        st.winv(st.w==0)=0;
        
        GhG=[];
        
    case 'ReversedCartesian'
        ParOrder = pe_indices(:,2);
        
        [tempy, tempz]=ndgrid(1:Ny*2,1:Nz);
        tempy=fftshift(tempy);
        tempz=fftshift(tempz);
        temp=sub2ind([Ny*2 Nz],tempy,tempz);
        
        st.pes=sub2ind([Ny*2 Nz],pe_indices(:,1)*2, ParOrder(:));
        st.pes=temp(st.pes);
        
        st.w=zeros(Ny*2*Nz,1);
        st.Ahmat=spalloc(Ny*2*Nz,numel(st.pes),numel(st.pes));
        for j=1:Ny*2*Nz
            st.w(j)=sum(st.pes==j);
            st.Ahmat(j,:)=(st.pes.'==j);
        end
        
        st.w(1)=st.w(1)-size(nav_data,1);
        st.winv = 1./st.w;
        st.winv(st.w==0)=0;
        
        GhG=[];
end