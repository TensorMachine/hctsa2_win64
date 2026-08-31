function HCTSA2_DiagnoseTiseanMex
% HCTSA2_DiagnoseTiseanMex   Show exactly how a TISEAN MEX output differs from
% the executable's, byte for byte.
%
%   HCTSA2_DiagnoseTiseanMex
%
% Run from the hctsa 2.0 root, after startup, with BOTH backends built. Prints
% a lot; copy the whole thing back.
%
% This exists because HCTSA2_TestTiseanMex reports only that a comparison
% failed. On a machine where the four Fortran programs disagree and the six C
% programs do not, the useful question is what the differing bytes actually are
% -- a line-ending difference, a formatting difference and a numeric difference
% look identical in a pass/fail count and have completely different causes.

if exist('BF_WriteTempFile','file') ~= 2
    error('HCTSA2_DiagnoseTiseanMex:noPath','Run startup first.');
end
fprintf('\n============ TISEAN MEX diagnosis ============\n');
fprintf('%s on %s\n',version,computer);
fprintf('tempdir : %s (%u chars)\n',tempdir,numel(tempdir));
t = tempname;
fprintf('tempname: %s (%u chars)\n',t,numel(t));
if numel(t) + 2 > 72
    fprintf(2,['  NOTE: over 72 characters. The Fortran programs declare their\n' ...
               '  filenames character*72 and truncate anything longer.\n']);
end

rng(1);
y = cumsum(randn(600,1)); y = (y-mean(y))/std(y);
fn = BF_WriteTempFile(y);
fprintf('data    : %s (%u chars)\n',fn,numel(fn));

% Short output paths, deliberately: this rules the 72-character limit in or out
% independently of everything else.
shortA = [fn '.a'];
shortB = [fn '.b'];

fprintf('\n--- 1. c2t, short output paths ---\n');
c2 = [fn '.c2'];
run1('d2',sprintf('-d1 -M1,3 -t10 -o %s %s',fn,fn),'exe');
if ~isfile(c2)
    fprintf(2,'  d2 did not produce %s -- cannot continue\n',c2); return
end
compareBackends('c2t',sprintf('-o %%s %s',c2),shortA,shortB);

fprintf('\n--- 2. c2t, tempname output paths (as the test harness uses) ---\n');
compareBackends('c2t',sprintf('-o %%s %s',c2),[tempname '_e'],[tempname '_m']);

fprintf('\n--- 3. c1, short output paths ---\n');
compareBackends('c1',sprintf('-d1 -m1 -M3 -t10 -n50 -o %%s %s',fn),shortA,shortB);

fprintf('\n--- 4. c1 twice on the SAME backend (is each backend self-consistent?) ---\n');
selfConsistency('c1',sprintf('-d1 -m1 -M3 -t10 -n50 -o %%s %s',fn),shortA,shortB,'exe');
selfConsistency('c1',sprintf('-d1 -m1 -M3 -t10 -n50 -o %%s %s',fn),shortA,shortB,'mex');

fprintf('\n--- 5. d2, short output paths (a C program, for contrast) ---\n');
compareBackends('d2',sprintf('-d1 -M1,3 -t10 -o %%s %s',fn),shortA,shortB,'.c2');

fprintf('\n--- 6. the strings BF_RunTisean returns (what the operations parse) ---\n');
% NL_TISEAN_c1 does not read c2d's file; it parses the string BF_RunTisean
% returns. On the exe backend that is captured stdout, on the MEX backend it is
% a temporary -o file read back. This compares the two directly.
run1('c1',sprintf('-d1 -m1 -M3 -t10 -n50 -o %s.c1 %s',fn,fn),'exe');
compareStrings('c2d',sprintf('-a2 %s.c1',fn));
if isfile([fn '.c2'])
    compareStrings('c2g',sprintf('%s.c2',fn));
    compareStrings('c2t',sprintf('%s.c2',fn));
end
compareStrings('nstat_z',sprintf('-# 4 -d1 -m3 %s',fn));

for e = {'','.c1','.c2','.d2','.h2','.stat','.a','.b','.a.c2','.b.c2','.tso'}
    f = [fn e{1}];
    if isfile(f), try delete(f); catch, end, end
end
fprintf('\n=============================================\n\n');
end

%-------------------------------------------------------------------------------
function compareStrings(prog,args)
setenv('HCTSA2_TISEAN_BACKEND','exe'); clear BF_RunTisean
[s1,a] = BF_RunTisean(sprintf('%s %s',prog,args));
setenv('HCTSA2_TISEAN_BACKEND','mex'); clear BF_RunTisean
[s2,b] = BF_RunTisean(sprintf('%s %s',prog,args));
fprintf('  %-8s status %d vs %d, %u vs %u chars',prog,s1,s2,numel(a),numel(b));
if isequal(a,b), fprintf('  -- IDENTICAL\n'); return; end
fprintf('  -- DIFFER\n');
na = numel(a); nb = numel(b); n = min(na,nb);
d = find(a(1:n) ~= b(1:n),1);
if isempty(d)
    fprintf(2,'    common prefix identical; lengths differ by %u\n',abs(na-nb));
    fprintf(2,'    exe tail: |%s|\n',escape(a(min(na,n+1):min(na,n+80))));
    fprintf(2,'    mex tail: |%s|\n',escape(b(min(nb,n+1):min(nb,n+80))));
else
    fprintf(2,'    first difference at character %u\n',d);
    lo = max(1,d-50);
    fprintf(2,'    exe: |%s|\n',escape(a(lo:min(na,d+50))));
    fprintf(2,'    mex: |%s|\n',escape(b(lo:min(nb,d+50))));
end
va = str2double(strsplit(strtrim(a))); vb = str2double(strsplit(strtrim(b)));
va = va(isfinite(va)); vb = vb(isfinite(vb));
fprintf(2,'    parsed %u vs %u numbers\n',numel(va),numel(vb));
if ~isempty(va) && ~isempty(vb)
    m = min(numel(va),numel(vb));
    fprintf(2,'    max relative difference over the common %u: %.3g\n', ...
            m,max(abs(va(1:m)-vb(1:m))./max(1e-30,abs(va(1:m)))));
end
end

function s = escape(t)
s = t;
s(s == char(10)) = '$';
s(s == char(13)) = '^';
s(s < 32) = '.';
end

%-------------------------------------------------------------------------------
function compareBackends(prog,argsFmt,outExe,outMex,ext)
if nargin < 5, ext = ''; end
run1(prog,sprintf(argsFmt,outExe),'exe');
run1(prog,sprintf(argsFmt,outMex),'mex');
a = [outExe ext]; b = [outMex ext];
fprintf('  exe -> %s\n',a);
fprintf('  mex -> %s\n',b);
describe(a,b);
tidy({a,b});
end

function selfConsistency(prog,argsFmt,o1,o2,backend)
run1(prog,sprintf(argsFmt,o1),backend);
run1(prog,sprintf(argsFmt,o2),backend);
fprintf('  %s backend, two runs: ',backend);
describe(o1,o2);
tidy({o1,o2});
end

%-------------------------------------------------------------------------------
function run1(prog,args,backend)
setenv('HCTSA2_TISEAN_BACKEND',backend);
clear BF_RunTisean
try
    BF_RunTisean(sprintf('%s %s',prog,args),true);
catch err
    fprintf(2,'    %s (%s) errored: %s\n',prog,backend,err.message);
end
end

%-------------------------------------------------------------------------------
function describe(a,b)
if ~isfile(a) || ~isfile(b)
    fprintf(2,'    MISSING: exists(a)=%d exists(b)=%d\n',isfile(a),isfile(b));
    return
end
x = readBytes(a); y = readBytes(b);
fprintf('    sizes %u vs %u',numel(x),numel(y));
if isequal(x,y), fprintf('  -- IDENTICAL\n'); return; end
fprintf('  -- DIFFER\n');

% Is it only line endings?
sx = x(x ~= 13); sy = y(y ~= 13);
if isequal(sx,sy)
    fprintf(2,'    the ONLY difference is carriage returns: %u in a, %u in b\n', ...
            sum(x==13),sum(y==13));
    return
end
fprintf('    CR count %u vs %u, LF count %u vs %u\n', ...
        sum(x==13),sum(y==13),sum(x==10),sum(y==10));

n = min(numel(x),numel(y));
d = find(x(1:n) ~= y(1:n),1);
if isempty(d)
    fprintf('    common prefix identical; one file is longer by %u bytes\n', ...
            abs(numel(x)-numel(y)));
    d = n;
end
fprintf('    first difference at byte %u\n',d);
lo = max(1,d-60); hi = min(n,d+60);
fprintf('    exe: |%s|\n',printable(x(lo:hi)));
fprintf('    mex: |%s|\n',printable(y(lo:min(numel(y),d+60))));

% Numeric comparison of the first column pair, if it parses.
va = parseNums(a); vb = parseNums(b);
if ~isempty(va) && ~isempty(vb)
    m = min(numel(va),numel(vb));
    fprintf('    parsed %u vs %u numbers',numel(va),numel(vb));
    if m > 0
        rel = abs(va(1:m)-vb(1:m)) ./ max(1e-30,abs(va(1:m)));
        fprintf('; max relative difference over the common %u: %.3g\n',m,max(rel));
    else
        fprintf('\n');
    end
end
end

function v = parseNums(f)
v = [];
try
    txt = fileread(f);
    c = textscan(txt,'%f','CommentStyle','#','Whitespace',' \t\r\n');
    v = c{1};
    v = v(isfinite(v));
catch
end
end

function b = readBytes(f)
fid = fopen(f,'r'); b = fread(fid,Inf,'*uint8'); fclose(fid); b = b(:)';
end

function s = printable(b)
s = char(b);
s(b < 32 | b > 126) = '.';
s = strrep(s,char(0),'.');
end

function tidy(files)
for i = 1:numel(files)
    for e = {'','.c1','.c2','.d2','.h2','.stat'}
        f = [files{i} e{1}];
        if isfile(f), try delete(f); catch, end, end
    end
end
end
