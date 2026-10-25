function ok = HCTSA2_AuditTiseanCleanup(beVocal)
% HCTSA2_AuditTiseanCleanup   Check every TISEAN temp file has a guarded exit.
%
%   ok = HCTSA2_AuditTiseanCleanup
%
% Reads the source of every operation that calls BF_WriteTempFile and checks,
% statically, that the temporary file is removed on EVERY way out. Returns true
% when they all pass. It is a scan of the source text for the pattern below, not a
% proof: it cannot see what a line does, only that the pattern is there.
%
%---WHY THIS IS A STATIC CHECK AND NOT A RUNTIME ONE
% A runtime test only covers the paths the test data happens to reach. Three
% separate attempts at this cleanup shipped with leaks precisely because each
% one enumerated the exits by hand and missed a kind:
%
%   attempt 1  onCleanup alone. Does not fire in a file that contains a NESTED
%              function -- the parent workspace outlives the return -- and six
%              of the ten operations of the time had one.
%   attempt 2  explicit calls before "return" and at the end of the body. Missed
%              all 43 error() calls.
%   attempt 3  a try/catch round the body. try/catch does not intercept RETURN,
%              and nine returns were written "out = NaN; return", on one line,
%              so a regex looking for return alone on a line found none of them.
%
% Every one of those was a whole CLASS of exit missed, not a typo, and each was
% found by the user on their machine rather than here. So the rule is checked
% mechanically now, over the real source, every time.
%
%---WHAT IT REQUIRES OF EACH OPERATION
%   1. the line after "filePath = BF_WriteTempFile(...)" opens a try block
%   2. that block ends with a catch that cleans up and rethrows
%   3. every return inside it cleans up first
%   4. the LAST statement of the block, the one before the catch, is the cleanup
%      (HCTSA2 2.2.1: this used to accept a cleanup on any line without a return,
%      which a cleanup inside a return branch also satisfied.)
%
% try/catch covers error() and anything a called function throws; the per-return
% calls cover what try/catch cannot see. Together they are exhaustive: there is
% no third way out of a MATLAB function.

if nargin < 1 || isempty(beVocal), beVocal = true; end

% Discover the callers rather than listing them, so an operation added later is
% audited without anyone remembering to add it here.
opDir = fileparts(which('NL_d2'));
listing = dir(fullfile(opDir,'*.m'));
ok = true; nChecked = 0;

if beVocal, fprintf('\n-- TISEAN temporary-file cleanup audit --\n'); end
for k = 1:numel(listing)
    f = fullfile(opDir,listing(k).name);
    src = strsplit(fileread(f),char(10));
    w = find(~cellfun(@isempty,regexp(src,'^\s*filePath\s*=\s*BF_WriteTempFile\(','once')),1);
    if isempty(w), continue; end
    nChecked = nChecked + 1;
    problems = auditOne(src,w);
    if isempty(problems)
        if beVocal, fprintf('  %-24s ok\n',listing(k).name); end
    else
        ok = false;
        fprintf(2,'  %-24s %s\n',listing(k).name,strjoin(problems,'; '));
    end
end
if nChecked == 0
    % Nothing found means the scan looked in the wrong place (startup not run,
    % so NL_d2 is not on the path), not that every caller is guarded.
    ok = false;
    fprintf(2,'  no TISEAN caller found in %s -- run startup first\n',opDir);
end
if beVocal
    fprintf('  %u file(s) checked, %s\n',nChecked, ...
            ternary(ok,'all guarded','PROBLEMS ABOVE'));
end
end

%-------------------------------------------------------------------------------
function problems = auditOne(src,w)
problems = {};

t = firstCode(src,w+1);
if isempty(t) || ~strcmp(strtrim(src{t}),'try')
    problems{end+1} = 'no try block opens after BF_WriteTempFile';
    return
end
c = find(~cellfun(@isempty,regexp(src,'^\s*catch\s+tiseanCleanupErr','once')),1);
if isempty(c) || c < t
    problems{end+1} = 'no matching "catch tiseanCleanupErr"';
    return
end
% Comments are dropped first: a comment that mentions BF_CleanTempFiles and
% rethrow must not pass for the code.
tail = strjoin(cellfun(@stripComment,src(c:min(c+4,numel(src))),'UniformOutput',false),' ');
if isempty(strfind(tail,'BF_CleanTempFiles')) %#ok<STREMP>
    problems{end+1} = 'catch does not clean up';
end
if isempty(strfind(tail,'rethrow')) %#ok<STREMP>
    problems{end+1} = 'catch does not rethrow';
end

for i = t+1:c-1
    code = stripComment(src{i});
    if isempty(strtrim(code)), continue; end
    hasClean  = ~isempty(strfind(code,'BF_CleanTempFiles(filePath)')); %#ok<STREMP>
    hasReturn = ~isempty(regexp(code,'\<return\>','once'));
    if hasReturn && ~hasClean
        % allow the cleanup to sit on the previous executable line
        prev = prevCode(src,i-1,t);
        if isempty(prev) || isempty(strfind(stripComment(src{prev}),'BF_CleanTempFiles(filePath)')) %#ok<STREMP>
            problems{end+1} = sprintf('line %u: return with no cleanup',i); %#ok<AGROW>
        end
    end
end
% The normal path: the last executable statement of the try body is the cleanup,
% at the nesting depth of the try itself, and is not a return.
last = prevCode(src,c-1,t+1);
lastCode = '';
if ~isempty(last), lastCode = stripComment(src{last}); end
if isempty(strfind(lastCode,'BF_CleanTempFiles(filePath)')) || ... %#ok<STREMP>
        ~isempty(regexp(lastCode,'\<return\>','once'))
    problems{end+1} = 'the last statement of the try block is not the cleanup';
end
end

%-------------------------------------------------------------------------------
function s = stripComment(line)
% Drop a trailing comment so "% ... return ..." never counts as a return.
q = false; s = line;
for i = 1:numel(line)
    ch = line(i);
    if ch == '''', q = ~q; end
    if ch == '%' && ~q, s = line(1:i-1); return; end
end
end

function i = firstCode(src,from)
i = [];
for k = from:numel(src)
    if ~isempty(strtrim(stripComment(src{k}))), i = k; return; end
end
end

function i = prevCode(src,from,floorIdx)
i = [];
for k = from:-1:max(floorIdx,1)
    if ~isempty(strtrim(stripComment(src{k}))), i = k; return; end
end
end

function v = ternary(c,a,b)
if c, v = a; else, v = b; end
end
