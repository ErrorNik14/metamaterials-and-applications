function make_output_folders(P)
%MAKE_OUTPUT_FOLDERS  Create results/, results/figures/, results/data/ if missing.
folders = {P.paths.results, P.paths.figs, P.paths.data, P.paths.anim};
for k = 1:numel(folders)
    if ~exist(folders{k}, 'dir')
        mkdir(folders{k});
    end
end
end
