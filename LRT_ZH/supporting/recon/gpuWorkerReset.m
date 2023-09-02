function gpuWorkerReset(varargin)
spmd
  try
    gpuDevice(mod(labindex-1,gpuDeviceCount)+1);
  catch %if fails, defaults to first gpu
  end
end
end