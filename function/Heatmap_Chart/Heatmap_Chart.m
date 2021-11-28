%%
% *This is an example of how to create a heatmap chart in MATLAB&#174;* .
% 
% You can open this example in the <https://www.mathworks.com/products/matlab/live-editor.html 
% Live Editor> with MATLAB version 2016a or higher.
%
% Read about the <http://www.mathworks.com/help/matlab/ref/heatmap.html |heatmap|> function in the MATLAB documentation. This function is available in R2017a or newer.
% For more examples, go to <http://www.mathworks.com/discovery/gallery.html MATLAB Plot Gallery>
%
% Copyright 2017-2018 The MathWorks, Inc.

% Check version
if verLessThan('matlab','9.2')
    error('heatmap is available in R2017a or newer.')
end

% Load ride data from Boston's bike sharing program
load CambridgeData cambridge

% Create a heatmap of DayOfWeek vs. AgeGroup, with color representing count
hm = heatmap(cambridge,'AgeGroup','DayOfWeek');

% Change the color to represent average Duration
hm.ColorVariable = 'Duration';
hm.ColorMethod = 'mean';
