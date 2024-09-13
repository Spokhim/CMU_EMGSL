function fig = plot_mesh_arm(vertices, faces, vertexNormals, vertexLabels, faceLabels, options)

arguments
    vertices (:,3) double
    faces (:,:) double
    vertexNormals (:,3) double
    vertexLabels (:,1) uint8
    faceLabels (:,1) uint8
    options.Alpha (3,1) double = [1;    % Bone
                                  0.8;  % Muscle
                                  0.3]; % Vessels
    options.Colors (3,3) double = [0.8, 0.8, 0.8;   % Gray for Bone
                                  1.0, 0.6, 0.6;   % Light red for Muscle
                                  0.6, 0.6, 1.0];  % Light blue for Vessel
    options.Labels (1,3) uint8 = uint8([1, 2, 3]);
end

% Create a figure for the 3D visualization
fig = figure('Color','w','Name','3D Arm Mesh Visualizer', ...
    'Units','inches','Position',[1.5 0.5 8 6]);
ax = axes(fig,'NextPlot','add','FontName','Tahoma','FontSize',14);

% Define colors for each segment (bone, muscle, vessel)
colors = options.Colors;
alpha = options.Alpha;

% Plot each segment separately
uniqueLabels = options.Labels;

for i = 1:length(uniqueLabels)
    % Get the label for the current segment
    label = uniqueLabels(i);
    
    % Find the faces that belong to this label
    facesIdx = faceLabels == label;
    verticesIdx = vertexLabels == label;
    firstVertexOffset = find(verticesIdx,1,'first')-1;
    
    % Plot the mesh for this segment
    patch(ax, ...
        'Vertices', vertices(verticesIdx, :), ...
        'Faces', faces(facesIdx, :)+1-firstVertexOffset, ...
        'FaceColor', colors(label, :), ...
        'FaceVertexCData', vertexNormals(verticesIdx,:), ...
        'EdgeColor', 'none', ...
        'FaceLighting', 'gouraud', ...
        'FaceAlpha', alpha(i));
end

% Set visualization properties
axis(ax, 'equal');
xlabel(ax,'X','FontName','Tahoma','FontSize',16);
ylabel(ax,'Y','FontName','Tahoma','FontSize',16);
zlabel(ax,'Z','FontName','Tahoma','FontSize',16);
title(ax, ...
    '3D Anatomy Visualization', ...
    'Bone, Muscle, and Vessel', ...
    'FontName','Tahoma','FontSize',20);
view(ax,3);
light(ax, 'Position', [0 0 10], 'Style', 'infinite');
material(ax, 'shiny');
camlight(ax);

% Add a legend for the segments
legend(ax,{'Bone', 'Muscle', 'Vessel'}, 'Location', 'northeastoutside');

end