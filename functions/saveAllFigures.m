function saveAllFigures(fileName)

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% savefig(FigList,[ fileName '.fig'])
savefig(FigList,sprintf('%s.fig', fileName));

end

