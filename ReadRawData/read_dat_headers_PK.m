function [headers,protocol]=read_dat_headers_PK(datfilename);
% function [headers,protocol]=read_dat_headers_PK(datfilename);
%
% function to read Siemens dat file headers for meas.dat (VB13)
%
% headers is output structure containing ascii protocol data
%     headers.Config
%     headers.Dicom
%     headers.MeasYaps
%     headers.Phoenix
%     headers.Spice

% protocol (optional output) is Matlab structure corresponding to MeasYaps
       
%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

fid=fopen(datfilename,'r');

hdrlen=fread(fid,1,'int32');

if isempty(hdrlen)
    headers=[];
    protocol=[];
    return
end

nbuffers=fread(fid,1,'int32');
for i=1:nbuffers
   tmp=1;
   name{i}=[];
   while tmp~=0;
       tmp=fread(fid,1,'char');
       name{i}=[name{i},tmp];
   end
   name{i}=char(name{i}(1:end-1));
   length(i)=fread(fid,1,'int32');
   tmp=char(fread(fid,length(i),'char')');
   eval(['headers.',name{i},'= tmp;']);   
end

if nargout==2
    protocol=MeasYaps2struct(headers);
end
fclose(fid);

% header description
% // |X X X X|X X X X|X X X X X....0|X X X X|X X X X X.....|X X X X X X....0|X X X X|X X X X X.....|XXXX..............|XXXX....
% // |4 bytes|4 bytes|   x Bytes    |4 Bytes|     x Bytes  |    x Bytes     |4 Bytes|  x Bytes     | padding          | data
% // |hdr len|buf nr.| name 1       |len 1  | prot 1       | name 2         |len 2  | prot 2       | 32byte aligned   |
% 
% Description:
% First 4 Bytes: Overall length of header (hdr len).
% This can be used to hop directly to the raw data.
% The next 4 bytes indicates the number of embedded data structures (buf nr.)
% Each data structure starts with a name (e.g. name 1):
% this is a NULL terminated string, then the next 4 bytes are
% indicating the length (e.g. len 1) of the data structure and
% the data structure (e.g. prot 1) itself.


function header=MeasYaps2struct(headers);
% function header=MeasYaps2struct(headers);
%
% function to covert MeasYaps from text to Matlab structure

%     ***************************************
%     *  Peter Kellman  (kellman@nih.gov)   *
%     *  Laboratory for Cardiac Energetics  *
%     *  NIH NHLBI                          *
%     ***************************************

MeasYaps=headers.MeasYaps;
lf = find(double(MeasYaps)==10); % find positions of line feeds
lf = [0,lf];

for line=1:length(lf)-1
    tline = MeasYaps(lf(line)+1:lf(line+1)-1);
    if ~ischar(tline); return; end
    index=findstr(tline,'=');
    
    if ~isempty(index)
        fieldname=deblank(tline(1:findstr(tline,'=')-1));
        fieldvalue=deblank(tline(findstr(tline,'=')+1:end));
        if max(isletter(fieldvalue))==0
            fieldvalue=str2num(fieldvalue);
        end
        fieldname=strrep(fieldname,'[','{');
        fieldname=strrep(fieldname,']','}');
        index2=findstr(fieldname,'{');
        index3=findstr(fieldname,'}');
 
%coding changes - vmpai - for vb12+ versions of IDEA. June 2005        
      
%         if ~isempty(index2)
%             fieldname = strrep(fieldname,fieldname(index2+1:index3-1),...
%                 num2str(str2num(fieldname(index2+1:index3-1))+1));
%         end

        eb = size(index2,2);
        ee = size(index3,2);
      
        strs = eb + 1;
        nums = eb;
        
        if ~isempty(index2)

            %split out the character portions of the statement
            str(1) = {fieldname(1:index2(1))};
            for i = 2:strs-1
                str(i) = {fieldname(index3(i-1):index2(i))};
            end
            str(strs) = {fieldname(index3(strs-1):end)};
            
			%pull out, increment and reconvert to string the numeric portion of the
			%statement
			for i=1:nums
                numstr(i) = {num2str(str2num(fieldname(index2(i)+1:index3(i)-1))+1)};
			end

			%combine the string and the numeric portions of the statement
			nstr = '';
            for i=1:eb
                nstr = strcat(nstr, str(i), numstr(i));
			end
            
           %convert from cell array to char array before pushing back in.
           fieldname = char(strcat(nstr,str(strs)));
%end coding change - vmpai - june 2005.           
        end
        eval(['header.',fieldname,'= fieldvalue;']);
    end
end

