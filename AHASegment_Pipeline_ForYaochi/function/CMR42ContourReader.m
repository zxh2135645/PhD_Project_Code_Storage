function con = CMR42ContourReader(cmr_file)
% This function reads CMR42 XML file and extracts contour points
%
%   con = CMR42ContourReader(cmr_file);
%
% Author: Avan Suinesiaputra - University of Auckland 2013

    try
        cmr42XML = xmlread(cmr_file);
    catch
        error('Cannot read %s file as an XML file.',cmr_file);
    end

    con.contours = [];

    % browse through items
    allItems = cmr42XML.getElementsByTagName('Hash:item');
    for k=0:allItems.getLength-1

        key = allItems.item(k).getAttribute('Hash:key');
        if( isempty(key) ), continue; end

        if( strcmpi(key,'Contours') )

            C = ReadContourNode(allItems.item(k));
            if( ~isempty(C) )
                con.contours = [con.contours; C];
            end

        elseif( strcmpi(key,'StudyUid') )

            con.studyIUID = char(allItems.item(k).getTextContent);

        end

    end

    fprintf(1,'Found %d contours\n',numel(con.contours));
end


% --- Reading points from one single contour node
function C = ReadContourNode(contour_node)
    
    C = [];
    
    % get number of contour types
    ntypes = str2num(contour_node.getAttribute('Hash:count'));
    if( ntypes < 1 )
        return; 
    end
    
    % this iuid should be the SOPInstanceUID
    C.iuid = char(contour_node.getParentNode.getAttribute('Hash:key'));
    C.ctype = {};
    C.pts = {};

    % get all the contours in this node
    subres = 1.0;
    for j=0:contour_node.getLength-1
        
        contourNode = contour_node.item(j);
        if( strcmpi(contourNode.getNodeName,'Hash:item') )
            
            % the name given as the type of the contour
            C.ctype = [C.ctype char(contourNode.getAttribute('Hash:key'))];
            
            % get all points
            for k=0:contourNode.getLength-1
                
                if( strcmpi(contourNode.item(k).getNodeName,'Hash:item') )
                    
                    % points
                    if( strcmpi(contourNode.item(k).getAttribute('Hash:key'),'Points') )
                    
                        pointNode = contourNode.item(k);
                        npoints = str2num(pointNode.getAttribute('List:count'));
                        pts = NaN * zeros(npoints,2);

                        xNodes = pointNode.getElementsByTagName('Point:x');
                        yNodes = pointNode.getElementsByTagName('Point:y');

                        if( xNodes.getLength ~= yNodes.getLength || xNodes.getLength ~= npoints )
                            error('Number of points mismatched.');
                        end

                        for p=0:npoints-1
                            pts(p+1,:) = [str2double(xNodes.item(p).getTextContent) str2double(yNodes.item(p).getTextContent)];
                        end


                        C.pts = [C.pts pts];
                        
                    % subpixel resolution
                    elseif( strcmpi(contourNode.item(k).getAttribute('Hash:key'),'SubpixelResolution') )
                        
                        subres = str2double(contourNode.item(k).getTextContent);
                        
                    end
                    
                end
                
            end
            
        end
    end
    
    if( numel(C.ctype) ~= ntypes )
        error('Number of contour types does not match with the number of items.');
    end
    
    % apply subres
    for i=1:numel(C.pts)
        C.pts{i} = C.pts{i} ./ subres;
    end
    
end