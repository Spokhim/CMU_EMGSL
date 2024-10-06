function fig = plot_annotated_template_slice(z, options)
% PLOT_ANNOTATED_TEMPLATE_SLICE Displays an annotated image of a specific anatomical slice.
% 
% This function reads a template image of a slice, loads associated landmarks
% from an Excel sheet, and optionally annotates the landmarks and image 
% on a specified or newly created axes.
%
% USAGE:
%   fig = plot_annotated_template_slice(z, options)
%
% INPUT:
%   z       - The slice index (positive integer) specifying which slice image and landmarks to load.
%
% OPTIONS (as name-value pairs):
%   FileExpression          - File path pattern to locate the slice images.
%                             Default: "C:/Data/Anatomy/Human Arm/Sections/R_Forearm_Section_%d.png"
%   LandMarksFile           - Excel file containing the landmarks for each slice.
%                             Default: "C:/Data/Anatomy/Human Arm/Sections/Landmarks.xlsx"
%   LandMarkSheetExpression - Sheet name pattern in the Excel file that contains the landmarks for the slice.
%                             Default: "R_Forearm_Section_%d"
%   LandmarksToAdd          - A string or array of strings specifying which landmarks to annotate. 
%                             Use "ALL" to annotate all landmarks.
%                             Default: "ALL"
%   AddLabels               - Logical flag to indicate if the landmarks should be labeled on the image.
%                             Default: true
%   AddImage                - Logical flag to indicate if the slice image should be displayed.
%                             Default: true
%   Axes                    - A handle to axes where the image and annotations will be plotted.
%                             If not provided, a new figure and axes are created.
%                             Default: []
%
% OUTPUT:
%   fig - A handle to the figure containing the annotated slice.
%
% DESCRIPTION:
%   This function plots a specific slice image from a given file path and 
%   annotates it with landmarks from an Excel file. It can either plot on 
%   an existing set of axes or create a new figure for the plot. If no axes 
%   are provided, a new figure is generated. Landmarks can be selectively 
%   annotated based on the 'LandmarksToAdd' option, and the function supports 
%   toggling of the image and labels with 'AddImage' and 'AddLabels'.
%
% EXAMPLE:
%   % Plot slice number 5 with default options
%   fig = plot_annotated_template_slice(5);
%
%   % Plot slice number 10 with custom landmarks and axes
%   ax = axes;
%   fig = plot_annotated_template_slice(10, 'LandmarksToAdd', ["A", "B"], 'Axes', ax);
%
% See also: IMREAD, READTABLE, AXES, TEXT

arguments
    z (1,1) {mustBePositive, mustBeInteger}
    options.FileExpression = "C:/Data/Anatomy/Human Arm/Sections/R_Forearm_Section_%d.png";
    options.LandMarksFile = "C:/Data/Anatomy/Human Arm/Sections/Landmarks.xlsx";
    options.LandMarkSheetExpression = "R_Forearm_Section_%d";
    options.LandmarksToAdd (1,:) string = "ALL";
    options.AddLabels (1,1) logical = true;
    options.AddImage (1,1) logical = true;
    options.Axes = [];
end

if isempty(options.Axes)
    fig = figure('Name','Annotated Template Arm Slice', ...
        'Color','w','Units','inches','Position',[2 2 7 5.5]);
    ax = axes(fig,'NextPlot','add','Color','none','YDir','reverse');
else
    ax = options.Axes;
    fig = ax.Parent;
end

I = imread(sprintf(options.FileExpression,z));
T = readtable(options.LandMarksFile, "Sheet", sprintf(options.LandMarkSheetExpression, z));

if options.AddImage
    imagesc(ax, I);
end
if options.AddLabels
    if strcmpi(options.LandmarksToAdd(1),"ALL")
        for ii = 1:size(T,1)
            text(ax, T.X(ii), T.Y(ii), T.Landmark{ii}, ...
                'FontName','Consolas','FontWeight','bold',...
                'FontSize',9, 'Color', 'w','HorizontalAlignment','center');
        end
    else
        for k = 1:numel(options.LandmarksToAdd)
            ii = ismember(T.Landmark,upper(options.LandmarksToAdd(k)));
            if nnz(ii)==1
                text(ax, T.X(ii), T.Y(ii), T.Landmark{ii}, ...
                    'FontName','Consolas','FontWeight','bold',...
                    'FontSize',9, 'Color', 'w','HorizontalAlignment','center');
            end
        end
    end
end


end