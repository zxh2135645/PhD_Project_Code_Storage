switch ScanType
  case 'T2star'
%     lr = 1e-4;
%     sms = 1e-6;  % Original: 1e-3; % How smooth should the wall_clock timecurve be?
%     smorders = 1;
%     circs = false;
    lr = 1e-4; %how low-rank?
    sms = [1e-6 1e-6 0];
    smorders = [1 1 1];
    circs = [true false false];
  case 'Cine'
    lr = 1e-4; %how low-rank?
    sms = 0;
  case 'IR'
    lr = 1e-4; %how low-rank?
    sms = [3e-6 1e-6];
    smorders = [1 1];
    circs = [true false];
  case 'T2IR'
    lr = 1e-4; %how low-rank?
    sms = [1e-6 1e-6 0];
    smorders = [1 1 1];
    circs = [true false false];
  case 'T2prep'
    lr = 1e-4;
    sms = [3e-6 1e-6 3e-6];
    smorders = [1 1 1];
    circs = [true false false];
  case 'SR'
    lr = 3e-5; %how low-rank?
    if max(Ridx)>1
      sms = [3e-6 1e-6 3e-6]; %how smooth? (cardiac, resp, then wall-clock dimension)
      smorders = [1 1 1]; %piecewise constant
      circs = [true false false];
    else
      sms = [3e-6 3e-6]; %how smooth? (cardiac, then wall-clock dimension)
      smorders = [1 1]; %piecewise constant
      circs = [true false];
    end
end

%Initial feature extraction
navdata_temp=reshape(reshape(nav_data,[],Ncoils)*mixer(:,1:ceil(Ncoils/3)),size(nav_data,1),[]);
[~,~,ts_proj]=svd(double(navdata_temp),'econ');
ts_proj=ts_proj(:,1:64);
navdata_temp=navdata_temp*ts_proj;