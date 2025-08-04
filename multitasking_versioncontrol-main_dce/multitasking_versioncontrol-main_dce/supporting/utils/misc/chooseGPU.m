for j=1:gpuDeviceCount
    gpudev = gpuDevice(j);
    am(j)  = gpudev.AvailableMemory;
end
[~,whichdev] = max(am);
ReconOptions.gpudev = gpuDevice(whichdev);