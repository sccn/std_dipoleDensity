% std_pop_dipoleDensity(): plots STUDY cluster dipoles and their density.
%
% See also: eegplugin_std_dipoleDensity() std_dipoleDensity() dipplot()        

% Author: Makoto Miyakoshi, JSPS/SCCN,INC,UCSD
% History
% 
% 10/02/2018 Makoto. Path to the Talailach tool added.
% 03/06/2017 Makoto. FWHM used. Default is changed to Zeynep and Luca's calculation to have 9.8 mm error in Gaussian.
% 02/24/2015 ver 0.22 by Makoto. (none) for the non-selected. BrainBlobBrower layout. 
% 01/16/2015 ver 0.20 by Makoto. Improved STUDY.design compatibility.
% 04/25/2014 ver 1.9 by Makoto. Change str2double -> str2num to avoid 'NaN' (Thanks Nikola Vukovic!)
% 05/22/2013 ver 1.8 by Makoto. Color red added.
% 03/29/2013 ver 1.7 by Makoto. Added cmin cmax. mir3dplot() is fixed accordingly.
% 02/19/2013 ver 1.6 by Makoto. STUDY = std_pop_dipoleDensity(). User inputs are now saved for the next use.
% 02/06/2013 ver 1.5 by Makoto. Color scale upper limit added.
% 01/28/2013 ver 1.4 by Makoto. Combined groups supported. Plot 1 - Plot 2 explained.
% 01/24/2013 ver 1.3 by Makoto. Group difference supported.
% 01/23/2013 ver 1.2 by Makoto. Single cluster supported.
% 12/03/2012 ver 1.1 by Makoto. Save figure added.
% 11/23/2012 ver 1.0 by Makoto. Created.

% Copyright (C) 2012, Makoto Miyakoshi JSPS/SCCN,INC,UCSD
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function STUDY = std_pop_dipoleDensity(STUDY, ALLEEG)

com = '';
if nargin < 2
    help std_pop_dipoleDensity;
    return;
end;

% Set path to 'dependency' folder.
currentFolderPath = which('std_pop_dipoleDensity');
addpathName       = [currentFolderPath(1:end-24) filesep 'dependency' filesep 'Talairach'];
addpath(addpathName)
disp('Path set: (this_plugin_root)/dependency/Talairach added.')

% create a letter string for clusters
clsString = '(none)|Parentcluster';
if length(STUDY.cluster)>1
    for n = 2:length(STUDY.cluster)
        clsString = [clsString '|'  STUDY.cluster(1,n).name];
    end
end

% group condition in the STUDY.design 
var2 = STUDY.design(STUDY.currentdesign).variable(1,2).value;
groupString = '(none)|all';
if ~isempty(var2)
    for n = 1:length(var2)
        if iscell(var2{1,n}) % if groups are combined in the selected group condition
            tmpCell = var2{1,n}; % work on copy
            tmpCell(2,:) = {' & '};
            tmpCell{2,end} = '';
            tmpString = [tmpCell{:}];
            clear tmpCell
        else
            tmpString = var2{1,n};
        end
        groupString = [groupString '|' tmpString];
    end
end   

% collect user input
try
   userInput = inputgui('title', 'std_pop_dipoleDensity()', 'geom', ...
       {{4 14 [0 0] [1 1]} {4 14 [1 0] [1 1]} {4 14 [2 0] [1 1]} {4 14 [3 0] [1 1]} ...
        {4 14 [0 1] [1 1]} {4 14 [1 1] [1 1]} {4 14 [2 1] [1 1]} {4 14 [3 1] [1 1]} ...
        {4 14 [0 2] [1 1]} {4 14 [1 2] [1 1]} {4 14 [2 2] [1 1]} {4 14 [3 2] [1 1]} ...
        {4 14 [0 3] [1 1]} {4 14 [1 3] [1 1]} {4 14 [2 3] [1 1]} {4 14 [3 3] [1 1]} ...
        {4 14 [0 4] [1 1]} {4 14 [1 4] [1 1]} {4 14 [2 4] [1 1]} {4 14 [3 4] [1 1]} ...
        {4 14 [0 5] [1 1]} {4 14 [1 5] [1 1]} {4 14 [2 5] [1 1]} {4 14 [3 5] [1 1]} ...
        {4 14 [1 7] [1 1]} {8 14 [4 7] [1 1]} ...
        {4 14 [1 8] [1 1]} {8 14 [4 8] [1 1]} ...
        {4 14 [1 9] [1 1]} {8 14 [4 9] [1 1]} ...
        {4 14 [1 11] [1 1]} {10 14 [4 11] [1 1]} {10 14 [5 11] [1 1]} {10 14 [6 11] [1 1]} {10 14 [7 11] [1 1]} {10 14 [8 11] [1 1]} ...
        {4 14 [1 14] [1 1]} {8 14 [4 14] [1 1]}}, ... 
    'uilist',...
       {{'style' 'text' 'string' ''} {'style' 'text' 'string' 'Cluster'} {'style' 'text' 'string' 'Group'} {'style' 'text' 'string' 'Color'} ...
        {'style' 'text' 'string' 'Plot 1: Cluster, Group, and Color'}    {'style' 'popupmenu' 'string' clsString 'tag' 'cluster' 'value' STUDY.dipplotWD{1,1}} {'style' 'popupmenu' 'string' groupString 'tag' 'group' 'value' STUDY.dipplotWD{1,2}} {'style' 'popupmenu' 'string' '(none)|White|Yellow|Fuchsia|Red|Silver|Gray|Olive|Purple|Maroon|Aqua|Lime|Teal|Green|Blue|Navy|Black' 'tag' 'color' 'value' STUDY.dipplotWD{1,3}} ...
        {'style' 'text' 'string' 'Plot 2: Cluster, Group, and Color'}    {'style' 'popupmenu' 'string' clsString 'tag' 'cluster' 'value' STUDY.dipplotWD{1,4}} {'style' 'popupmenu' 'string' groupString 'tag' 'group' 'value' STUDY.dipplotWD{1,5}} {'style' 'popupmenu' 'string' '(none)|White|Yellow|Fuchsia|Red|Silver|Gray|Olive|Purple|Maroon|Aqua|Lime|Teal|Green|Blue|Navy|Black' 'tag' 'color' 'value' STUDY.dipplotWD{1,6}} ...
        {'style' 'text' 'string' 'Plot 3: Cluster, Group, and Color'}    {'style' 'popupmenu' 'string' clsString 'tag' 'cluster' 'value' STUDY.dipplotWD{1,7}} {'style' 'popupmenu' 'string' groupString 'tag' 'group' 'value' STUDY.dipplotWD{1,8}} {'style' 'popupmenu' 'string' '(none)|White|Yellow|Fuchsia|Red|Silver|Gray|Olive|Purple|Maroon|Aqua|Lime|Teal|Green|Blue|Navy|Black' 'tag' 'color' 'value' STUDY.dipplotWD{1,9}} ...
        {'style' 'text' 'string' 'Plot 4: Cluster, Group, and Color'}    {'style' 'popupmenu' 'string' clsString 'tag' 'cluster' 'value' STUDY.dipplotWD{1,10}} {'style' 'popupmenu' 'string' groupString 'tag' 'group' 'value' STUDY.dipplotWD{1,11}} {'style' 'popupmenu' 'string' '(none)|White|Yellow|Fuchsia|Red|Silver|Gray|Olive|Purple|Maroon|Aqua|Lime|Teal|Green|Blue|Navy|Black' 'tag' 'color' 'value' STUDY.dipplotWD{1,12}} ...
        {'style' 'text' 'string' 'Plot 5: Cluster, Group, and Color'}    {'style' 'popupmenu' 'string' clsString 'tag' 'cluster' 'value' STUDY.dipplotWD{1,13}} {'style' 'popupmenu' 'string' groupString 'tag' 'group' 'value' STUDY.dipplotWD{1,14}} {'style' 'popupmenu' 'string' '(none)|White|Yellow|Fuchsia|Red|Silver|Gray|Olive|Purple|Maroon|Aqua|Lime|Teal|Green|Blue|Navy|Black' 'tag' 'color' 'value' STUDY.dipplotWD{1,15}} ...
        {'style' 'text' 'string' 'Slice orientation'}                    {'style' 'popupmenu' 'string' 'Axial|Sagittal|Coronal' 'tag' 'slice' 'value' STUDY.dipplotWD{1,16}} ...
        {'style' 'text' 'string' 'Smoothing Gaussian kernel FWHM [mm]'}  {'style' 'edit' 'string' STUDY.dipplotWD{1,17}} ...
        {'style' 'text' 'string' 'Color range [min max] blank=auto'}     {'style' 'edit' 'string' STUDY.dipplotWD{1,18}} ...
        {'style' 'text' 'string' 'If Plot 1 - Plot 2'} {'style' 'popupmenu' 'string' groupString 'value' STUDY.dipplotWD{1,19}} {'style' 'text' 'string' 'minus'} {'style' 'popupmenu' 'string' groupString 'value' STUDY.dipplotWD{1,20}} {'style' 'text' 'string' 'thresh. [%]'} {'style' 'edit' 'string' STUDY.dipplotWD{1,21}} ...
        {'style' 'text' 'string' 'Save figures to the current folder'}   {'style' 'checkbox', 'value', STUDY.dipplotWD{1,22}}});
catch
    userInput = inputgui('title', 'std_pop_dipoleDensity()', 'geom', ...
       {{4 14 [0 0] [1 1]} {4 14 [1 0] [1 1]} {4 14 [2 0] [1 1]} {4 14 [3 0] [1 1]} ...
        {4 14 [0 1] [1 1]} {4 14 [1 1] [1 1]} {4 14 [2 1] [1 1]} {4 14 [3 1] [1 1]} ...
        {4 14 [0 2] [1 1]} {4 14 [1 2] [1 1]} {4 14 [2 2] [1 1]} {4 14 [3 2] [1 1]} ...
        {4 14 [0 3] [1 1]} {4 14 [1 3] [1 1]} {4 14 [2 3] [1 1]} {4 14 [3 3] [1 1]} ...
        {4 14 [0 4] [1 1]} {4 14 [1 4] [1 1]} {4 14 [2 4] [1 1]} {4 14 [3 4] [1 1]} ...
        {4 14 [0 5] [1 1]} {4 14 [1 5] [1 1]} {4 14 [2 5] [1 1]} {4 14 [3 5] [1 1]} ...
        {4 14 [1 7] [1 1]} {8 14 [4 7] [1 1]} ...
        {4 14 [1 8] [1 1]} {8 14 [4 8] [1 1]} ...
        {4 14 [1 9] [1 1]} {8 14 [4 9] [1 1]} ...
        {4 14 [1 11] [1 1]} {10 14 [4 11] [1 1]} {10 14 [5 11] [1 1]} {10 14 [6 11] [1 1]} {10 14 [7 11] [1 1]} {10 14 [8 11] [1 1]} ...
        {4 14 [1 14] [1 1]} {8 14 [4 14] [1 1]}}, ... 
    'uilist',...
       {{'style' 'text' 'string' ''} {'style' 'text' 'string' 'Cluster'} {'style' 'text' 'string' 'Group'} {'style' 'text' 'string' 'Color'} ...
        {'style' 'text' 'string' 'Plot 1: Cluster, Group, and Color'}    {'style' 'popupmenu' 'string' clsString 'tag' 'cluster' 'value' 1} {'style' 'popupmenu' 'string' groupString 'tag' 'group' 'value' 1} {'style' 'popupmenu' 'string' '(none)|White|Yellow|Fuchsia|Red|Silver|Gray|Olive|Purple|Maroon|Aqua|Lime|Teal|Green|Blue|Navy|Black' 'tag' 'color' 'value' 1} ...
        {'style' 'text' 'string' 'Plot 2: Cluster, Group, and Color'}    {'style' 'popupmenu' 'string' clsString 'tag' 'cluster' 'value' 1} {'style' 'popupmenu' 'string' groupString 'tag' 'group' 'value' 1} {'style' 'popupmenu' 'string' '(none)|White|Yellow|Fuchsia|Red|Silver|Gray|Olive|Purple|Maroon|Aqua|Lime|Teal|Green|Blue|Navy|Black' 'tag' 'color' 'value' 1} ...
        {'style' 'text' 'string' 'Plot 3: Cluster, Group, and Color'}    {'style' 'popupmenu' 'string' clsString 'tag' 'cluster' 'value' 1} {'style' 'popupmenu' 'string' groupString 'tag' 'group' 'value' 1} {'style' 'popupmenu' 'string' '(none)|White|Yellow|Fuchsia|Red|Silver|Gray|Olive|Purple|Maroon|Aqua|Lime|Teal|Green|Blue|Navy|Black' 'tag' 'color' 'value' 1} ...
        {'style' 'text' 'string' 'Plot 4: Cluster, Group, and Color'}    {'style' 'popupmenu' 'string' clsString 'tag' 'cluster' 'value' 1} {'style' 'popupmenu' 'string' groupString 'tag' 'group' 'value' 1} {'style' 'popupmenu' 'string' '(none)|White|Yellow|Fuchsia|Red|Silver|Gray|Olive|Purple|Maroon|Aqua|Lime|Teal|Green|Blue|Navy|Black' 'tag' 'color' 'value' 1} ...
        {'style' 'text' 'string' 'Plot 5: Cluster, Group, and Color'}    {'style' 'popupmenu' 'string' clsString 'tag' 'cluster' 'value' 1} {'style' 'popupmenu' 'string' groupString 'tag' 'group' 'value' 1} {'style' 'popupmenu' 'string' '(none)|White|Yellow|Fuchsia|Red|Silver|Gray|Olive|Purple|Maroon|Aqua|Lime|Teal|Green|Blue|Navy|Black' 'tag' 'color' 'value' 1} ...
        {'style' 'text' 'string' 'Slice orientation'}                    {'style' 'popupmenu' 'string' 'Axial|Sagittal|Coronal' 'tag' 'slice' 'value' 1} ...
        {'style' 'text' 'string' 'Smoothing Gaussian kernel FWHM [mm]'}  {'style' 'edit' 'string' '14.2'} ...
        {'style' 'text' 'string' 'Color range [min max] blank=auto'}     {'style' 'edit' 'string' ''} ...
        {'style' 'text' 'string' 'If Plot 1 - Plot 2'}                   {'style' 'popupmenu' 'string' groupString 'value' 1} {'style' 'text' 'string' 'minus'} {'style' 'popupmenu' 'string' groupString 'value' 1} {'style' 'text' 'string' 'thresh. [%]'} {'style' 'edit' 'string' '5'} ...
        {'style' 'text' 'string' 'Save figures to the current folder'}   {'style' 'checkbox', 'value', 0}});
end

% canceled
if isempty(userInput)
    return
end

% store userInput to STUDY
STUDY.dipplotWD = userInput;

% create colorIndex
colorIndex = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];

% 16 colors names officially supported by W3C specification for HTML
colors{1,1}  = [1 1 1];            % White
colors{2,1}  = [1 1 0];            % Yellow
colors{3,1}  = [1 0 1];            % Fuchsia
colors{4,1}  = [1 0 0];            % Red
colors{5,1}  = [0.75  0.75  0.75]; % Silver
colors{6,1}  = [0.5 0.5 0.5];      % Gray
colors{7,1}  = [0.5 0.5 0];        % Olive
colors{8,1}  = [0.5 0 0.5];        % Purple
colors{9,1}  = [0.5 0 0];          % Maroon
colors{10,1} = [0 1 1];            % Aqua
colors{11,1} = [0 1 0];            % Lime
colors{12,1} = [0 0.5 0.5];        % Teal
colors{13,1} = [0 0.5 0];          % Green
colors{14,1} = [0 0 1];            % Blue
colors{15,1} = [0 0 0.5];          % Navy
colors{16,1} = [0 0 0];            % Black

colorName{1,1}  = 'White';
colorName{2,1}  = 'Yellow';
colorName{3,1}  = 'Fuchsia';
colorName{4,1}  = 'Red';
colorName{5,1}  = 'Silver';
colorName{6,1}  = 'Gray';
colorName{7,1}  = 'Olive';
colorName{8,1}  = 'Purple';
colorName{9,1}  = 'Maroon';
colorName{10,1} = 'Aqua';
colorName{11,1} = 'Lime';
colorName{12,1} = 'Teal';
colorName{13,1} = 'Green';
colorName{14,1} = 'Blue';
colorName{15,1} = 'Navy';
colorName{16,1} = 'Black';

plotParams{1,1}.cluster = userInput{1,1};
plotParams{1,1}.group   = userInput{1,2};
if userInput{1,3}>1
    plotParams{1,1}.color     = colors{colorIndex(userInput{1,3}-1)};
    plotParams{1,1}.colorName = colorName{colorIndex(userInput{1,3}-1)};
end

plotParams{1,2}.cluster = userInput{1,4};
plotParams{1,2}.group   = userInput{1,5};
if userInput{1,6}>1
    plotParams{1,2}.color     = colors{colorIndex(userInput{1,6}-1)};
    plotParams{1,2}.colorName = colorName{colorIndex(userInput{1,6}-1)};
end

plotParams{1,3}.cluster = userInput{1,7};
plotParams{1,3}.group   = userInput{1,8};
if userInput{1,9}>1
    plotParams{1,3}.color     = colors{colorIndex(userInput{1,9}-1)};
    plotParams{1,3}.colorName = colorName{colorIndex(userInput{1,9}-1)};
end

plotParams{1,4}.cluster = userInput{1,10};
plotParams{1,4}.group   = userInput{1,11};
if userInput{1,12}>1
    plotParams{1,4}.color     = colors{colorIndex(userInput{1,12}-1)};
    plotParams{1,4}.colorName = colorName{colorIndex(userInput{1,12}-1)};
end

plotParams{1,5}.cluster = userInput{1,13};
plotParams{1,5}.group   = userInput{1,14};
if userInput{1,15}>1
    plotParams{1,5}.color     = colors{colorIndex(userInput{1,15}-1)};
    plotParams{1,5}.colorName = colorName{colorIndex(userInput{1,15}-1)};
end

plotParams{1,6}  = userInput{1,16};             % slice orientation

% FWHM = 2.355*sigma See https://en.wikipedia.org/wiki/Full_width_at_half_maximum
plotParams{1,7}  = str2num(userInput{1,17})/2.355; %#ok<ST2NM> % Smoothing [mm] converted from Gauss sigma to FWHM

plotParams{1,8}  = str2num(userInput{1,18}); %#ok<ST2NM> % color scale upper limit [mm]
plotParams{1,9}  = userInput{1,19};             % group subtractor
plotParams{1,10} = userInput{1,20};             % group subtracted
plotParams{1,11} = str2num(userInput{1,21}); %#ok<ST2NM> % threshold [%]
plotParams{1,12} = userInput{1,22};             % save figure or not

% pass them to the main function
std_dipoleDensity(STUDY, ALLEEG, plotParams)