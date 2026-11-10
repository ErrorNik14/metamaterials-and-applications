function save_figure(fh, name, P)
%SAVE_FIGURE  Export a figure as a 300-dpi PNG (for slides) and optionally .fig.
%   Files go to results/figures/<name>.png and <name>.fig
png_path = fullfile(P.paths.figs, [name '.png']);
exportgraphics(fh, png_path, 'Resolution', P.plot.dpi, 'BackgroundColor', 'white');
if P.plot.save_fig
    savefig(fh, fullfile(P.paths.figs, [name '.fig']));
end
if P.plot.save_vector
    exportgraphics(fh, fullfile(P.paths.figs, [name '.pdf']), ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
end
fprintf('  saved %s\n', png_path);
end
