function report = HCTSA2_BuildBinaries(stageDir)
% HCTSA2_BuildBinaries   Compile every MEX component of hctsa 2.0.
%
%   HCTSA2_BuildBinaries
%   report = HCTSA2_BuildBinaries('mex_win64')
%
% Run this from the hctsa 2.0 root ON A MACHINE THAT HAS A C/C++ COMPILER.
% It compiles each component independently, never aborting the whole run
% because one target failed, and reports how many library features each
% success or failure is worth.
%
% If STAGEDIR is given, every binary produced is also copied there with its
% relative path preserved. That folder can then be carried to the offline
% machine and unpacked over the tree -- which is the intended workflow when the
% target machine has no compiler.
%
%---WHY THIS EXISTS
% 1,299 of the 7,702 library features (16.9%) need compiled code. The other
% 83.1% are pure MATLAB and work as soon as the paths are set. Nothing here is
% required to get a useful installation; each component simply unlocks its own
% slice.
%
%---ON WINDOWS
% MATLAB needs a supported compiler. The free option is MinGW-w64, normally
% installed through the Add-On Explorer, which needs internet. So on an offline
% target either (a) build on a connected machine and carry the staged binaries
% across, or (b) pre-install MinGW-w64 before going offline. Run `mex -setup`
% to confirm before starting.

if nargin < 1, stageDir = ''; end

root = pwd;
if exist(fullfile(root,'Toolboxes'),'dir') ~= 7
    error('HCTSA2_BuildBinaries:wrongDir', ...
        'Run this from the hctsa 2.0 root (the folder containing Toolboxes/).');
end

fprintf('\n=============== hctsa 2.0 binary build ===============\n');
fprintf('MATLAB %s on %s\n',version,computer);

% Compiler check up front -- failing here explains every later failure at once
haveC = true;
try
    cc = mex.getCompilerConfigurations('C','Selected');
    if isempty(cc)
        haveC = false;
    else
        fprintf('C compiler   : %s\n',cc(1).Name);
    end
catch
    haveC = false;
end
if ~haveC
    fprintf(2,'No C compiler is configured. Run "mex -setup" first.\n');
    fprintf(2,'Continuing anyway so the report shows what would be built.\n');
end
try
    cxx = mex.getCompilerConfigurations('C++','Selected');
    if ~isempty(cxx), fprintf('C++ compiler : %s\n',cxx(1).Name); end
catch
end
fprintf('\n');

% target list: {label, folder, files, featureCount, kind, note}
%
% kind is 'gates' (the features do not run without this MEX) or 'speed' (a
% working MATLAB fallback exists and the MEX only makes it faster). Everything
% is built either way -- a free speed-up is worth having -- but the summary
% keeps them separate so a 'speed' failure is not mistaken for lost features.
%
% solve_chol is the ONLY speed-only target in the package. Every other C source
% with a same-named .m (MS_complexitybs, MS_nearest, MS_shannon, nn_prepare)
% has a documentation-only stub with no executable code, so the MEX is required.
T = {
 'Max_Little/fastdfa',      fullfile('Toolboxes','Max_Little','fastdfa'),        {'ML_fastdfa_core.c'},                          10, 'gates', ''
 'Max_Little/rpde',         fullfile('Toolboxes','Max_Little','rpde'),           {'ML_close_ret.c'},                              0, 'gates', 'shares the Max_Little slice'
 'Max_Little/steps_bumps',  fullfile('Toolboxes','Max_Little','steps_bumps_toolkit'), {'ML_kvsteps_core.cpp'},                    0, 'gates', 'C++; shares the Max_Little slice'
 'Michael_Small',           fullfile('Toolboxes','Michael_Small'),               {'MS_complexitybs.c','MS_nearest.c','MS_shannon.c'}, 33, 'gates', ''
 'Physionet/sampen',        fullfile('Toolboxes','Physionet'),                   {'sampen_mex.c'},                               49, 'gates', 'gates EN_SampEn / EN_Randomize'
 'catch22',                 fullfile('Toolboxes','catch22','wrap_Matlab'),       {},                                            240, 'gates', 'uses its own mexAll'
 'gpml/solve_chol',         fullfile('Toolboxes','gpml','util'),                 {},                                             72, 'speed', 'Cholesky solve; MATLAB fallback exists, MEX is faster'
};

report = struct('name',{},'status',{},'features',{},'kind',{},'detail',{});
for i = 1:size(T,1)
    name = T{i,1}; rel = T{i,2}; files = T{i,3}; nfeat = T{i,4};
    kind = T{i,5}; note = T{i,6};
    fprintf('%-26s ',name);
    d = fullfile(root,rel);
    if exist(d,'dir') ~= 7
        fprintf('SKIP   (folder missing)\n');
        report(end+1) = mk(name,'missing',nfeat,'folder not present'); %#ok<AGROW>
        continue
    end
    ok = true; detail = '';
    cd(d);
    try
        if ~isempty(files)
            for k = 1:numel(files)
                if exist(fullfile(d,files{k}),'file') ~= 2
                    ok = false; detail = sprintf('%s not found',files{k}); break
                end
                mex(files{k});
            end
        else
            switch name
                case 'catch22'
                    mexAll;
                case 'gpml/solve_chol'
                    make;                   % compiles solve_chol.c
            end
        end
    catch err
        ok = false; detail = firstLine(err.message);
    end
    cd(root);
    if strcmp(kind,'speed'), tag = 'speed-up'; else, tag = 'features'; end
    if ok
        fprintf('OK     (%4d %-8s)  %s\n',nfeat,tag,note);
        report(end+1) = mk(name,'ok',nfeat,kind,note); %#ok<AGROW>
    else
        fprintf(2,'FAILED (%4d %-8s)  %s\n',nfeat,tag,detail);
        report(end+1) = mk(name,'failed',nfeat,kind,detail); %#ok<AGROW>
    end
end

%-------------------------------------------------------------------------------
%% Summary
%-------------------------------------------------------------------------------
isGate  = strcmp({report.kind},'gates');
isOK    = strcmp({report.status},'ok');
okF     = sum([report(isGate &  isOK).features]);
badF    = sum([report(isGate & ~isOK).features]);
spdOK   = sum([report(~isGate &  isOK).features]);
spdBad  = sum([report(~isGate & ~isOK).features]);
fprintf('\n---------------------------------------------------\n');
fprintf('features unlocked   : %5d\n',okF);
fprintf('features still gated: %5d\n',badF);
fprintf('pure MATLAB         : %5d (always available)\n',6403);
if spdOK > 0
    fprintf('accelerated         : %5d features now use MEX instead of the\n',spdOK);
    fprintf('                            MATLAB fallback (same results, faster)\n');
elseif spdBad > 0
    fprintf('not accelerated     : %5d features fall back to MATLAB (still correct)\n',spdBad);
end
fprintf('---------------------------------------------------\n');
fprintf('TISEAN (276 features) is NOT built here -- it is a set of standalone\n');
fprintf('executables, not MEX. See INSTALL.md section 4.\n');

%-------------------------------------------------------------------------------
%% Stage the binaries for transport
%-------------------------------------------------------------------------------
if ~isempty(stageDir)
    if exist(stageDir,'dir') ~= 7, mkdir(stageDir); end
    ext = ['*.' mexext];
    n = 0;
    L = dir(fullfile(root,'Toolboxes','**',ext));
    for i = 1:numel(L)
        relPath = strrep(L(i).folder,[root filesep],'');
        dest = fullfile(stageDir,relPath);
        if exist(dest,'dir') ~= 7, mkdir(dest); end
        copyfile(fullfile(L(i).folder,L(i).name),fullfile(dest,L(i).name));
        n = n + 1;
    end
    fprintf('\nStaged %u %s files into %s\n',n,mexext,stageDir);
    fprintf('Copy that folder onto the offline machine and merge it over the\n');
    fprintf('hctsa2 tree, preserving paths.\n');
end
fprintf('\n');
end

%-------------------------------------------------------------------------------
function s = mk(name,status,features,kind,detail)
s = struct('name',name,'status',status,'features',features,'kind',kind,'detail',detail);
end

function s = firstLine(msg)
s = strtrim(regexprep(msg,'\n.*','','once'));
if numel(s) > 70, s = [s(1:67) '...']; end
end
