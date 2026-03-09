%% Batch MEG preprocessing for all subjects under MEG_path
% Parallel version: parallelize across subjects, keep runs serial
%
% Pipeline:
% 1) Convert CTF .ds -> SPM M/EEG continuous
% 2) Mark artefacts on continuous data (MARK only, no rejection yet)
% 3) Epoch from marked continuous files
% 4) Manual baseline correction using [-3300 -2300] ms
% 5) Merge baseline-corrected epoched runs
%
% Notes:
% - Automatic baseline during epoching is disabled via S2.bc = 0
% - Artefacts are marked now and can be rejected later if needed
% - A log file is written into each subject output folder
% - This version is parfor-safe
% - Completed subjects are skipped using PROCESSING_COMPLETE.txt

clc; clear;

%% Initialise SPM
spm('Defaults','eeg');
spm_jobman('initcfg');

% disable popup
set(0,'DefaultFigureVisible','off')
%% User settings
MEG_path    = '/Volumes/ritd-ag-project-rd0270-eholm03/MEG';
output_root = fullfile(pwd, 'output');

prestim_sec   = 3.3;
poststim_sec  = 3.2;
baseline_ms   = [-3300 -2300];
EV_TYPE       = 'UPPT001_up';

% Artefact marking settings
artefact_threshold      = 5;
artefact_excwin_ms      = 100;
artefact_badchanthresh  = 1;

% Parallel workers
nWorkers = 4;   % adjust based on your Mac RAM/CPU

if ~isfolder(output_root)
    mkdir(output_root);
end

%% Load trial definition once
trialdef_file = fullfile(MEG_path, 'trialdef.mat');
if ~isfile(trialdef_file)
    error('trialdef.mat not found: %s', trialdef_file);
end

T  = load(trialdef_file);
td = T.trialdef;

code_num = arrayfun(@(x) x.eventvalue, td);
lab_str  = arrayfun(@(x) char(x.conditionlabel), td, 'UniformOutput', false);
lab_map  = containers.Map(num2cell(code_num), lab_str);

%% Find subject folders
subj_dirs = dir(MEG_path);
subj_dirs = subj_dirs([subj_dirs.isdir]);
subj_dirs = subj_dirs(~ismember({subj_dirs.name}, {'.','..'}));

% Keep numeric subject folders only, e.g. 201, 202, ...
is_num = cellfun(@(x) all(isstrprop(x, 'digit')), {subj_dirs.name});
subj_dirs = subj_dirs(is_num);

fprintf('Found %d subject folders under %s\n', numel(subj_dirs), MEG_path);

%% Skip subjects already completed
todo_mask = false(numel(subj_dirs), 1);

for s = 1:numel(subj_dirs)
    subj_no = subj_dirs(s).name;
    out_dir = fullfile(output_root, subj_no);
    done_file = fullfile(out_dir, 'PROCESSING_COMPLETE.txt');

    if ~isfile(done_file)
        todo_mask(s) = true;
    end
end

subj_dirs = subj_dirs(todo_mask);

fprintf('Subjects remaining to process: %d\n', numel(subj_dirs));
fprintf('Starting parallel pool with %d workers...\n', nWorkers);

if isempty(subj_dirs)
    fprintf('Nothing to do. All detected subjects are already completed.\n');
    return;
end

%% Start parallel pool
p = gcp('nocreate');
if isempty(p)
    parpool('local', nWorkers);
elseif p.NumWorkers ~= nWorkers
    delete(p);
    parpool('local', nWorkers);
end

%% Process all subjects in parallel
parfor s = 1:numel(subj_dirs)

    subj_no  = subj_dirs(s).name;
    subj_dir = fullfile(MEG_path, subj_no);
    out_dir  = fullfile(output_root, subj_no);

    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    done_file = fullfile(out_dir, 'PROCESSING_COMPLETE.txt');
    log_file  = fullfile(out_dir, sprintf('processing_subject_%s.txt', subj_no));

    % Worker ID for clearer logs
    task = getCurrentTask();
    if isempty(task)
        wid = 0;
    else
        wid = task.ID;
    end

    % Safety: skip if somehow already completed
    if isfile(done_file)
        logmsg(log_file, 'Worker %d | Subject %s already completed. Skipping.', wid, subj_no);
        continue;
    end

    if exist(log_file, 'file')
        delete(log_file);
    end

    try
        logmsg(log_file, '====================================================');
        logmsg(log_file, 'Worker %d | Processing subject %s (%d/%d)', wid, subj_no, s, numel(subj_dirs));
        logmsg(log_file, '====================================================');
        logmsg(log_file, 'Start time: %s', datestr(now));

        %% Handling files
        runs = dir(fullfile(subj_dir, '*.ds'));
        runs = runs([runs.isdir]);

        [~, idx] = sort({runs.name});
        runs = runs(idx);
        nRuns = numel(runs);

        if nRuns == 0
            logmsg(log_file, 'Worker %d | WARNING: Subject %s has no .ds runs. Skipping.', wid, subj_no);
            logmsg(log_file, 'Finished subject %s at %s', subj_no, datestr(now));
            continue;
        end

        logmsg(log_file, 'Worker %d | Found %d runs for subject %s', wid, nRuns, subj_no);

        %% 1) Convert runs
        D_files = cell(nRuns,1);

        t_convert = tic;
        for i = 1:nRuns
            ds_path   = fullfile(runs(i).folder, runs(i).name);
            run_label = sprintf('sub-%s_run-%02d', subj_no, i);

            logmsg(log_file, 'Worker %d | (%02d/%02d) Converting %s', wid, i, nRuns, run_label);

            S = [];
            S.dataset = ds_path;
            S.mode    = 'continuous';
            S.outfile = fullfile(out_dir, [run_label '.mat']);

            D = spm_eeg_convert(S);
            D_files{i} = fullfile(D.path, D.fname);
        end
        logmsg(log_file, 'Worker %d | All runs converted for subject %s.', wid, subj_no);
        logmsg(log_file, 'Worker %d | Conversion time: %.2f sec', wid, toc(t_convert));

        %% 2) Artefact marking on continuous data (MARK only)
        A_files = cell(nRuns,1);

        t_art = tic;
        for i = 1:nRuns

            if isempty(D_files{i}) || ~isfile(D_files{i})
                logmsg(log_file, 'Worker %d | WARNING: Subject %s Run %02d converted file missing; skipping artefact marking.', wid, subj_no, i);
                continue;
            end

            D = spm_eeg_load(D_files{i});
            logmsg(log_file, 'Worker %d | Subject %s Run %02d: artefact marking on %s', wid, subj_no, i, D.fname);

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

            try
                Da = spm_eeg_artefact(S);
                A_files{i} = fullfile(Da.path, Da.fname);
                logmsg(log_file, 'Worker %d | Subject %s Run %02d: artefact-marked file = %s', wid, subj_no, i, Da.fname);
            catch ME_art
                A_files{i} = '';
                logmsg(log_file, 'Worker %d | WARNING: Subject %s Run %02d artefact marking failed: %s', wid, subj_no, i, ME_art.message);
                logmsg(log_file, 'Worker %d | Falling back to converted file for epoching on this run.', wid);
            end
        end
        logmsg(log_file, 'Worker %d | Done: artefact marking finished for subject %s.', wid, subj_no);
        logmsg(log_file, 'Worker %d | Artefact marking time: %.2f sec', wid, toc(t_art));

        %% 3) Epoching from artefact-marked continuous files
        E_files = cell(nRuns,1);

        t_epoch = tic;
        for i = 1:nRuns

            if ~isempty(A_files{i}) && isfile(A_files{i})
                run_mat = A_files{i};
            elseif ~isempty(D_files{i}) && isfile(D_files{i})
                run_mat = D_files{i};
            else
                logmsg(log_file, 'Worker %d | WARNING: Subject %s Run %02d no input file found for epoching.', wid, subj_no, i);
                continue;
            end

            D = spm_eeg_load(run_mat);

            if ~strcmpi(D.type, 'continuous')
                logmsg(log_file, 'Worker %d | WARNING: Subject %s Run %02d is not continuous (D.type=%s). Skipping.', wid, subj_no, i, D.type);
                continue;
            end

            fs    = D.fsample;
            preS  = round(prestim_sec  * fs);
            postS = round(poststim_sec * fs);

            ev = events(D);
            if iscell(ev)
                ev = [ev{:}];
            end

            % parfor-safe explicit initialisation
            esamp = [];
            etime = [];

            if istable(ev)
                etype = string(ev.type);
                hasSample = ismember('sample', ev.Properties.VariableNames);
                hasTime   = ismember('time',   ev.Properties.VariableNames);

                if hasSample
                    esamp = ev.sample;
                end
                if hasTime
                    etime = ev.time;
                end

                vals_raw = ev.value;
                if ~iscell(vals_raw)
                    vals_raw = num2cell(vals_raw);
                end
            else
                etype = string({ev.type});
                hasSample = isfield(ev, 'sample');
                hasTime   = isfield(ev, 'time');

                if hasSample
                    esamp = [ev.sample];
                end
                if hasTime
                    etime = [ev.time];
                end

                vals_raw = {ev.value};
            end

            idxType = (etype == string(EV_TYPE));
            vals_raw_type = vals_raw(idxType);

            vnum = nan(numel(vals_raw_type),1);
            for k = 1:numel(vals_raw_type)
                x = vals_raw_type{k};

                if isempty(x)
                    vnum(k) = NaN;
                elseif isnumeric(x) && isscalar(x)
                    vnum(k) = double(x);
                elseif ischar(x) || isstring(x)
                    vnum(k) = str2double(strtrim(string(x)));
                elseif iscell(x) && numel(x)==1
                    y = x{1};
                    if isnumeric(y) && isscalar(y)
                        vnum(k) = double(y);
                    elseif ischar(y) || isstring(y)
                        vnum(k) = str2double(strtrim(string(y)));
                    else
                        vnum(k) = NaN;
                    end
                else
                    vnum(k) = NaN;
                end
            end

            kept_codes = unique(vnum(~isnan(vnum)));
            logmsg(log_file, 'Worker %d | Subject %s Run %02d event codes found: %s', wid, subj_no, i, mat2str(kept_codes'));

            idxCode = ismember(vnum, code_num);
            idxKeep = find(idxType);
            idxKeep = idxKeep(idxCode);

            if isempty(idxKeep)
                logmsg(log_file, 'Worker %d | WARNING: Subject %s Run %02d no %s events with codes of interest found.', wid, subj_no, i, EV_TYPE);
                continue;
            end

            if ~isempty(esamp) && numel(esamp) == numel(etype)
                samp = esamp(idxKeep);
            elseif ~isempty(etime) && numel(etime) == numel(etype)
                samp = round(etime(idxKeep) * fs) + 1;
            else
                logmsg(log_file, 'Worker %d | WARNING: Subject %s Run %02d events() has neither usable sample nor time field.', wid, subj_no, i);
                continue;
            end

            codes = vnum(idxCode);

            [samp, ord] = sort(samp(:));
            codes = codes(ord);

            beg = samp - preS;
            fin = samp + postS;

            valid = (beg >= 1) & (fin <= D.nsamples);

            nDropped = sum(~valid);
            beg   = beg(valid);
            fin   = fin(valid);
            codes = codes(valid);

            if isempty(beg)
                logmsg(log_file, 'Worker %d | WARNING: Subject %s Run %02d all trials invalid after boundary check.', wid, subj_no, i);
                continue;
            end

            offset = -preS * ones(size(beg));
            trl = [beg fin offset];

            condlabels = cell(size(codes));
            for kk = 1:numel(codes)
                condlabels{kk} = lab_map(codes(kk));
            end

            if nDropped > 0
                logmsg(log_file, 'Worker %d | Subject %s Run %02d: trials defined = %d (dropped %d boundary trials)', ...
                    wid, subj_no, i, size(trl,1), nDropped);
            else
                logmsg(log_file, 'Worker %d | Subject %s Run %02d: trials defined = %d', wid, subj_no, i, size(trl,1));
            end

            S2 = struct();
            S2.D = D;
            S2.trl = trl;
            S2.conditionlabels = condlabels;
            S2.prefix = 'e';
            S2.bc = 0;

            Dep = spm_eeg_epochs(S2);

            logmsg(log_file, 'Worker %d | Subject %s Run %02d: epoched ntrials = %d (%s)', wid, subj_no, i, Dep.ntrials, Dep.fname);
            E_files{i} = fullfile(Dep.path, Dep.fname);
        end
        logmsg(log_file, 'Worker %d | Done: each run epoched for subject %s.', wid, subj_no);
        logmsg(log_file, 'Worker %d | Epoching time: %.2f sec', wid, toc(t_epoch));

        %% 4) Manual baseline correction
        B_files = cell(nRuns,1);

        t_bc = tic;
        for i = 1:nRuns

            if isempty(E_files{i}) || ~isfile(E_files{i})
                logmsg(log_file, 'Worker %d | WARNING: Subject %s Run %02d epoched file missing.', wid, subj_no, i);
                continue;
            end

            D = spm_eeg_load(E_files{i});

            S = [];
            S.D = D;
            S.timewin = baseline_ms;

            Db = spm_eeg_bc(S);

            logmsg(log_file, 'Worker %d | Subject %s Run %02d: baseline corrected (%s)', wid, subj_no, i, Db.fname);
            B_files{i} = fullfile(Db.path, Db.fname);
        end
        logmsg(log_file, 'Worker %d | Done: baseline correction finished for subject %s.', wid, subj_no);
        logmsg(log_file, 'Worker %d | Baseline correction time: %.2f sec', wid, toc(t_bc));

        %% 4b) Optional baseline check
        idx_check = find(~cellfun(@isempty, E_files) & ~cellfun(@isempty, B_files), 1, 'first');
        if ~isempty(idx_check)
            De = spm_eeg_load(E_files{idx_check});
            Db = spm_eeg_load(B_files{idx_check});

            ch = indchantype(De, 'MEG');
            if isempty(ch)
                ch = [indchantype(De,'MEGMAG') indchantype(De,'MEGPLANAR')];
            end

            if ~isempty(ch)
                ch = ch(1);
                ix = (De.time >= baseline_ms(1)/1000) & (De.time <= baseline_ms(2)/1000);
                xe = squeeze(De(ch,:,1));
                xb = squeeze(Db(ch,:,1));

                logmsg(log_file, 'Worker %d | Baseline check on subject %s run %02d', wid, subj_no, idx_check);
                logmsg(log_file, 'Worker %d | Channel: %d (%s)', wid, ch, De.chanlabels{ch});
                logmsg(log_file, 'Worker %d | Before manual BC, baseline mean = %.12f', wid, mean(xe(ix)));
                logmsg(log_file, 'Worker %d | After  manual BC, baseline mean = %.12f', wid, mean(xb(ix)));
            end
        end

        %% 5) Merge baseline-corrected epoched runs
        valid_B = B_files(~cellfun(@isempty, B_files));

        if numel(valid_B) < 2
            logmsg(log_file, 'Worker %d | WARNING: Subject %s needs at least 2 baseline-corrected files to merge. Found %d.', ...
                wid, subj_no, numel(valid_B));
        else
            t_merge = tic;

            cwd0 = pwd;
            cd(out_dir);

            S = [];
            S.D = char(valid_B{:});
            S.recode = 'same';
            S.prefix = 'c';

            D_merged = spm_eeg_merge(S);

            cd(cwd0);

            merged_file = fullfile(D_merged.path, D_merged.fname);
            logmsg(log_file, 'Worker %d | Subject %s merged file saved at: %s', wid, subj_no, merged_file);

            Dm = spm_eeg_load(merged_file);
            logmsg(log_file, 'Worker %d | Subject %s merged dataset: %d trials, %d channels, %d samples per trial', ...
                wid, subj_no, Dm.ntrials, Dm.nchannels, Dm.nsamples);

            logmsg(log_file, 'Worker %d | Merge time: %.2f sec', wid, toc(t_merge));

            % completion marker: only write after successful end-to-end processing
            fid = fopen(done_file, 'w');
            if fid ~= -1
                fprintf(fid, 'Subject %s completed at %s\n', subj_no, datestr(now));
                fclose(fid);
                logmsg(log_file, 'Worker %d | Completion marker written: %s', wid, done_file);
            else
                logmsg(log_file, 'Worker %d | WARNING: Could not write completion marker: %s', wid, done_file);
            end
        end

        logmsg(log_file, 'Worker %d | Finished subject %s at %s', wid, subj_no, datestr(now));

    catch ME
        logmsg(log_file, 'Worker %d | ERROR while processing subject %s', wid, subj_no);
        logmsg(log_file, 'Worker %d | Message: %s', wid, ME.message);
        for ee = 1:numel(ME.stack)
            logmsg(log_file, 'Worker %d |   In %s at line %d', wid, ME.stack(ee).name, ME.stack(ee).line);
        end
    end
end

fprintf('\nAll subjects processed.\n');

%% Local logging helper
function logmsg(log_file, fmt, varargin)
    msg = sprintf(fmt, varargin{:});
    timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    line = sprintf('[%s] %s\n', timestamp, msg);

    % print to terminal
    fprintf('%s', line);

    % append to subject log file
    fid = fopen(log_file, 'a');
    if fid ~= -1
        fprintf(fid, '%s', line);
        fclose(fid);
    else
        fprintf(2, 'Could not open log file: %s\n', log_file);
    end
end