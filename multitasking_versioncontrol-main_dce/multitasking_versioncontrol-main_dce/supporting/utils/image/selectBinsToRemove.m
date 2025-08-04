function selectedBins = selectBinsToRemove(Nbins,img)

Items = {'Bin 1','Bin 2','Bin 3','Bin 4','Bin 5','Bin 6','Bin 7','Bin 8','Bin 9','Bin 10'};
Items = Items(1:Nbins);

if nargin < 2 || size(squeeze(img),3) < 2
    fig = uifigure('Name','Select Bins to Remove','Position',[100 100 240 320]);
    g = uigridlayout(fig,[Nbins+2 1]);
    temp = cell(1,Nbins+2);
    temp(:) = {'fit'};
    temp(Nbins+1) = {20};
    g.RowHeight = temp;
    g.RowSpacing = 8;
    g.ColumnWidth = {'1x'};
    
    cbx = cell(Nbins,1);
    for n = 1:Nbins
        cbx{n} = uicheckbox(g,"Text",Items(n));
        cbx{n}.Layout.Row = n;
    end
    
    % Confirm button
    ButtonH = uibutton(g,"state",'Text','Confirm Selection');
    ButtonH.Layout.Row = Nbins+2;
    waitfor(ButtonH,'Value');
    
    selectedBins = zeros(1,Nbins);
    for n = 1:Nbins
        selectedBins(n) = cbx{n}.Value;
    end
    selectedBins = find(selectedBins);
    
    close(fig);
else
    fig = uifigure('Name','Select Resp Bins to Remove');
    g = uigridlayout(fig,[Nbins+4 2]);
    temp = cell(1,Nbins+4);
    temp(:) = {'1x'};
    g.RowHeight = temp;
    g.ColumnWidth = {'3x','1x'};
    ax = uipanel(g);
    ax.Layout.Row = [1 Nbins+4];
    ax.Layout.Column = 1;
    s = sliceViewer(squeeze(img),'Parent',ax);  
    
    titleLabel = uilabel(g,'Text','Select bins to remove');
    titleLabel.Layout.Row = 1;

    cbx = cell(Nbins,1);
    for n = 1:Nbins
        cbx{n} = uicheckbox(g,"Text",Items(n));
        cbx{n}.Layout.Row = n+2;
        cbx{n}.Layout.Column = 2;
    end    

    ButtonH = uibutton(g,"state",'Text','Confirm Selection');
    ButtonH.Layout.Row = Nbins+4;
    waitfor(ButtonH,'Value');
    selectedBins = zeros(1,Nbins);
    for n = 1:Nbins
        selectedBins(n) = cbx{n}.Value;
    end
    selectedBins = find(selectedBins);
    close(fig); 
end
