%% Simple single-subject MEG preprocessing with SPM
clc; clear;

%% Initialise SPM
spm('Defaults','eeg');
spm_jobman('initcfg');

%% User settings
MEG_path    = '[your_data_directory]';
output_root = '[workspace_directory]';

subj_no = '201';   % change subject here, choosed 201 for example

prestim_sec = 3.3;
poststim_sec = 3.2;
baseline_ms = [-3300 -2300];
EV_TYPE = 'UPPT001_up';

artefact_threshold = 5;
artefact_excwin_ms = 100;
artefact_badchanthresh = 1;

%% Paths
subj_dir = fullfile(MEG_path, subj_no);
out_dir  = fullfile(output_root, subj_no);

if ~isfolder(out_dir)
    mkdir(out_dir);
end

%% Load trial definition
T  = load(fullfile(MEG_path, 'trialdef.mat'));
td = T.trialdef;

code_num = arrayfun(@(x) x.eventvalue, td);
lab_str  = arrayfun(@(x) char(x.conditionlabel), td, 'UniformOutput', false);
lab_map  = containers.Map(num2cell(code_num), lab_str);

%% Find runs
runs = dir(fullfile(subj_dir, '*.ds'));
runs = runs([runs.isdir]);

[~, idx] = sort({runs.name});
runs = runs(idx);
nRuns = numel(runs);

%% 1) Convert runs
D_files = cell(nRuns,1);

for i = 1:nRuns
    ds_path   = fullfile(runs(i).folder, runs(i).name);
    run_label = sprintf('sub-%s_run-%02d', subj_no, i);

    S = [];
    S.dataset = ds_path;
    S.mode    = 'continuous';
    S.outfile = fullfile(out_dir, [run_label '.mat']);

    D = spm_eeg_convert(S);
    D_files{i} = fullfile(D.path, D.fname);
end

%% 2) Artefact marking on continuous data
A_files = cell(nRuns,1);

for i = 1:nRuns
    D = spm_eeg_load(D_files{i});

    S = [];
    S.D      = D;
    S.mode   = 'mark';
    S.prefix = 'a';

    S.methods = struct([]);
    S.methods(1).channels = {'MEG'};
    S.methods(1).fun = 'zscorediff';
    S.methods(1).settings.threshold = artefact_threshold;
    S.methods(1).settings.excwin = artefact_excwin_ms;
    S.methods(1).settings.badchanthresh = artefact_badchanthresh;

    Da = spm_eeg_artefact(S);
    A_files{i} = fullfile(Da.path, Da.fname);
end

%% 3) Epoching
E_files = cell(nRuns,1);

for i = 1:nRuns
    D = spm_eeg_load(A_files{i});

    fs    = D.fsample;
    preS  = round(prestim_sec  * fs);
    postS = round(poststim_sec * fs);

    ev = events(D);
    if iscell(ev)
        ev = [ev{:}];
    end

    if istable(ev)
        etype = string(ev.type);
        esamp = ev.sample;
        vals_raw = ev.value;
        if ~iscell(vals_raw)
            vals_raw = num2cell(vals_raw);
        end
    else
        etype = string({ev.type});
        esamp = [ev.sample];
        vals_raw = {ev.value};
    end

    idxType = (etype == string(EV_TYPE));
    vals_raw = vals_raw(idxType);

    vnum = nan(numel(vals_raw),1);
    for k = 1:numel(vals_raw)
        x = vals_raw{k};

        if isnumeric(x) && isscalar(x)
            vnum(k) = double(x);
        elseif ischar(x) || isstring(x)
            vnum(k) = str2double(strtrim(string(x)));
        elseif iscell(x) && numel(x)==1
            y = x{1};
            if isnumeric(y) && isscalar(y)
                vnum(k) = double(y);
            else
                vnum(k) = str2double(strtrim(string(y)));
            end
        end
    end

    idxCode = ismember(vnum, code_num);
    idxKeep = find(idxType);
    idxKeep = idxKeep(idxCode);

    samp  = esamp(idxKeep);
    codes = vnum(idxCode);

    [samp, ord] = sort(samp(:));
    codes = codes(ord);

    beg = samp - preS;
    fin = samp + postS;

    offset = -preS * ones(size(beg));
    trl = [beg fin offset];

    condlabels = cell(size(codes));
    for kk = 1:numel(codes)
        condlabels{kk} = lab_map(codes(kk));
    end

    S2 = [];
    S2.D = D;
    S2.trl = trl;
    S2.conditionlabels = condlabels;
    S2.prefix = 'e';
    S2.bc = 0;

    Dep = spm_eeg_epochs(S2);
    E_files{i} = fullfile(Dep.path, Dep.fname);
end

%% 4) Baseline correction
B_files = cell(nRuns,1);

for i = 1:nRuns
    D = spm_eeg_load(E_files{i});

    S = [];
    S.D = D;
    S.timewin = baseline_ms;

    Db = spm_eeg_bc(S);
    B_files{i} = fullfile(Db.path, Db.fname);
end

%% 5) Merge baseline-corrected runs
cwd0 = pwd;
cd(out_dir);

S = [];
S.D = char(B_files{:});
S.recode = 'same';
S.prefix = 'c';

D_merged = spm_eeg_merge(S);

cd(cwd0);
