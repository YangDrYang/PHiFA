function h = plotTarget(cmptarget)
colors = tubscolors();
h = figure('Name','Target Segment Object Display');
hold on;

segcolors(1,:) = colors.GreenLight;
segcolors(2,:) = colors.VioletLight;
segcolors(3,:) = colors.BlueLight;
segcolors(4,:) = colors.Orange;
segcolors(5,:) = colors.tubsRed;
segcolors(6,:) = colors.BlueDark;
segcolors(7,:) = colors.BlueDark;

for j = 1:length(cmptarget.segments)
    for i = 1:length(cmptarget.segments(j).facets)
        if isa(cmptarget.segments(j).facets(i),'clRectangle')
            vtri = [cmptarget.segments(j).offset+cmptarget.segments(j).facets(i).base...
                cmptarget.segments(j).offset+cmptarget.segments(j).facets(i).base+cmptarget.segments(j).facets(i).edge1...
                cmptarget.segments(j).offset+cmptarget.segments(j).facets(i).base+cmptarget.segments(j).facets(i).edge1+cmptarget.segments(j).facets(i).edge2...
                cmptarget.segments(j).offset+cmptarget.segments(j).facets(i).base+cmptarget.segments(j).facets(i).edge2].';
            f = [1 2 3 4];
            patch('Faces',f,'Vertices',vtri,'FaceColor',       segcolors(j,:), ...
                 'EdgeColor',       'black',        ...
                 'FaceLighting',    'flat',     ...
                 'FaceAlpha', 0.2);
        else
            vtri = [cmptarget.segments(j).offset+cmptarget.segments(j).facets(i).base...
                cmptarget.segments(j).offset+cmptarget.segments(j).facets(i).base+cmptarget.segments(j).facets(i).edge1...
                cmptarget.segments(j).offset+cmptarget.segments(j).facets(i).base+cmptarget.segments(j).facets(i).edge2].';
            f = [1 2 3];
            patch('Faces',f,'Vertices',vtri,'FaceColor',       segcolors(j,:), ...
                 'EdgeColor',       'black',        ...
                 'FaceLighting',    'flat',     ...
                 'FaceAlpha', 0.2);
        end
    end
end
view([30 20]);
xlabel('x [m]','Fontsize',14);
ylabel('y [m]','Fontsize',14);
zlabel('z [m]','Fontsize',14);
% set(gca,'xticklabel',[])
% set(gca,'yticklabel',[])
% set(gca,'zticklabel',[])
axis equal;
grid;
end

