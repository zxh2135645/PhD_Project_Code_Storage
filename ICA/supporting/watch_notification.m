function watch_notification(mainpath,message)
if isunix
  sep = '/';
else
  sep = '\';
end
fid=fopen(strcat(mainpath,sep,'watch',sep,message),'w');
fprintf(fid,'%s',datestr(now));
fclose(fid);
return