addpath('/Users/hojungcho/Downloads/spm');

%% Finding .ds directory
data_dir = "../ds003082/sub-0001";

dsDirs = dir(fullfile(data_dir, '**', '*.ds'));
dsDirs = dsDirs([dsDirs.isdir]);

fifFiles = dir(fullfile(data_dir, '**', '*_meg.fif'));

fprintf('========= MEG DATASETS =========\n\n');

fprintf('CTF .ds folders found: %d\n', numel(dsDirs));
for i = 1:numel(dsDirs)
    fprintf('  %s\n', fullfile(dsDirs(i).folder, dsDirs(i).name));
end
fprintf('\n');

fprintf('FIF MEG files found: %d\n', numel(fifFiles));
for i = 1:numel(fifFiles)
    fprintf('  %s\n', fullfile(fifFiles(i).folder, fifFiles(i).name));
end
fprintf('\n');

fprintf('========= DONE =========\n');

%% 