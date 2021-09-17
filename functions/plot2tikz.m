function plot2tikz(name, width, varargin)

if ~exist('tex', 'dir')
    mkdir('tex');
end

print2 = sprintf('tex/%s.tex',name);
width = sprintf('%.1f\\textwidth',width);
if nargin>2
    extraaxop = sprintf('%s', varargin{1});
    matlab2tikz(print2, 'width', width, 'showWarnings', false, 'showInfo', ...
        false, 'extraAxisOptions', extraaxop, 'floatFormat', '%.5f');
else
    matlab2tikz(print2, 'width', width, 'showWarnings', false, 'showInfo', ...
        false, 'floatFormat', '%.5f');
end

end

