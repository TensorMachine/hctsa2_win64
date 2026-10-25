function HCTSA2_TiseanProbe(y)
% HCTSA2_TiseanProbe   Capture what the TISEAN executables actually print.
%
%   HCTSA2_TiseanProbe            % uses a short synthetic series
%   HCTSA2_TiseanProbe(y)         % uses your own
%
% Shows stdout and stderr SEPARATELY, which is the distinction that matters:
% TISEAN writes data to stdout and its banner to stderr, and MATLAB's system()
% merges them in an order that differs between platforms.

% (HCTSA2 2.2.1: rng(1) below re-seeds the global stream; leave it as the caller
% had it.)
rngState = rng;
restoreRng = onCleanup(@() rng(rngState)); %#ok<NASGU>
if nargin < 1 || isempty(y)
    rng(1); y = cumsum(randn(1000,1)); y = (y-mean(y))/std(y);
end
y = y(:);

fprintf('\n=========== TISEAN output probe ===========\n');
fprintf('MATLAB %s on %s\n',version,computer);

fprintf('\n--- which executables are visible ---\n');
% (HCTSA2 2.2.1: asks whether the file is there, in the folders BF_RunTisean
% looks in, then on the PATH. It used to run `<name> -h` and take exit code 127 or
% 9009 to mean "not found", which is the POSIX shell's code, not necessarily
% cmd.exe's, and several programs exit non-zero on -h anyway.)
cmds = {'d2','nstat_z','false_nearest','boxcount','lyap_r', ...
        'lyap_spec','nrlazy','c1','c2d','c2g','c2t','ripser'};
for i = 1:numel(cmds)
    [where,how] = locate(cmds{i});
    if isempty(where)
        verdict = 'NOT FOUND';
    else
        verdict = sprintf('found (%s): %s',how,where);
    end
    fprintf('  %-16s %s\n',cmds{i},verdict);
end

% An hc###### file from BF_WriteTempFile, removed with whatever the programs write
% beside it (.c2, .d2, .h2, .stat) on every way out. The same '%15.6f' format.
fn = BF_WriteTempFile(y,'%15.6f');
cleanProbe = onCleanup(@() BF_CleanTempFiles(fn)); %#ok<NASGU>
fprintf('\ndata file: %s  (%u chars)\n',fn,numel(fn));

if ispc, nul = 'nul'; else, nul = '/dev/null'; end

cmd = sprintf('nstat_z -# 4 -d1 -m3 "%s"',fn);
fprintf('\n--- nstat_z ---\n');
runBoth(cmd,nul);

cmd = sprintf('d2 -d1 -M1,5 -t50 "%s"',fn);
fprintf('\n--- d2 ---\n');
runBoth(cmd,nul);
fprintf('  files produced: ');
for e = {'.c2','.d2','.h2'}
    fprintf('%s=%d ',e{1},exist([fn e{1}],'file')==2);
end
fprintf('\n');

if exist([fn '.c2'],'file') == 2
    fprintf('\n--- c2g (reads d2''s .c2 output) ---\n');
    runBoth(sprintf('c2g "%s.c2"',fn),nul);
end

fprintf('\n=========== end of probe ===========\n\n');
end

%-------------------------------------------------------------------------------
function [where,how] = locate(prog)
% Where PROG is, by file: the folders startup and BF_RunTisean use (the same
% rule as HCTSA2_CheckInstall), then anywhere on the PATH. Empty if nowhere.
root = fileparts(fileparts(mfilename('fullpath'))); % (HCTSA2 3.0.0: this file is in Tests/)
if ispc
    ext = '.exe';
    dirs = {fullfile(root,'Toolboxes','Tisean_bin')};
else
    ext = '';
    tiseanDir = fullfile(root,'Toolboxes','Tisean_3.0.1','bin');
    if exist(tiseanDir,'dir') ~= 7, tiseanDir = fullfile(getenv('HOME'),'bin'); end
    dirs = {tiseanDir, fullfile(root,'Toolboxes','ripser','bin')};
end
where = ''; how = '';
for i = 1:numel(dirs)
    f = fullfile(dirs{i},[prog ext]);
    if isfile(f), where = f; how = 'install folder'; return; end
end
if ispc, probe = ['where ' prog ext ' 2>nul']; else, probe = ['command -v ' prog ' 2>/dev/null']; end
[st,w] = system(probe);
if st == 0
    where = strtrim(strtok(w,char(10))); how = 'on the PATH only';
end
end

%-------------------------------------------------------------------------------
function runBoth(cmd,nul)
[~,merged] = system(cmd);                         % what hctsa used to see
[~,clean]  = system([cmd ' 2>' nul]);             % stdout only
fprintf('  MERGED  (stdout+stderr): %u chars, %u lines\n', ...
    numel(merged),numel(strsplit(merged,char(10))));
show(merged,8);
fprintf('  STDOUT ONLY            : %u chars, %u lines\n', ...
    numel(clean),numel(strsplit(clean,char(10))));
show(clean,8);
end

function show(res,n)
lines = strsplit(res,char(10));
for i = 1:min(n,numel(lines))
    L = lines{i};
    if any(L == char(13))
        tag = '   <-- has CR';
    else
        tag = '';
    end
    fprintf('    %2u |%s|%s\n',i,strrep(L,char(13),''),tag);
end
if numel(lines) > n
    fprintf('    ... %u more lines\n',numel(lines)-n);
end
end
