function available_memory=getGPUmem(varargin)
[~,d]=unix(sprintf('nvidia-smi | grep MiB | head -%0d | awk ''{print $9, $11}''',gpuDeviceCount));
d=textscan(d,'%d%s%d%s');
available_memory=d{3}-d{1};
end
