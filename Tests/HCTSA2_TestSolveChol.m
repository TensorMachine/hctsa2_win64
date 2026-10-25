function HCTSA2_TestSolveChol
% HCTSA2_TestSolveChol   Is the gpml solve_chol MEX returning correct answers?
%
%   HCTSA2_TestSolveChol
%
% solve_chol(L,B) must return A\B where A = L'*L. gpml's infExact uses it for
% alpha and for the covariance derivatives, so a wrong result there does not
% announce itself -- it corrupts the gradient, the optimiser wanders into a
% region where the covariance matrix is no longer positive definite, and the
% failure finally surfaces as
%
%   Inference method failed [Error using chol / Matrix must be positive definite]
%
% a long way from the cause. This checks the MEX directly against MATLAB's own
% A\B, for the right-hand sides gpml uses (a matrix, a vector, eye(n)). If only
% the pure-MATLAB solve_chol.m is there, there is no MEX to test and it says so.
%
% Run after startup.

fprintf('\n=========== solve_chol check ===========\n');
w = which('solve_chol');
fprintf('resolving to : %s\n',w);
switch exist('solve_chol','file')
    case 3, fprintf('type         : MEX (compiled)\n');
    case 2
        % (HCTSA2 2.2.1: with only the .m there is no MEX in play, and comparing
        % the .m with A\B and then saying "the MEX agrees" would be false.)
        if isfile(fullfile(fileparts(w),['solve_chol.' mexext]))
            fprintf('type         : MEX built beside the .m\n');
        else
            fprintf('type         : .m (MATLAB fallback -- the MEX is not built)\n');
            fprintf('Nothing to test: there is no solve_chol MEX. The 148 gpml features run on\n');
            fprintf('the MATLAB implementation. Build it with HCTSA2_BuildBinaries.\n');
            fprintf('=======================================\n\n');
            return
        end
    otherwise, fprintf(2,'solve_chol not found on the path.\n'); return
end

% (HCTSA2 2.2.1: rng(1) re-seeds the global stream; leave it as the caller had it.)
rngState = rng;
restoreRng = onCleanup(@() rng(rngState)); %#ok<NASGU>
rng(1);
sizes = [5 20 100 300];
fprintf('\n%6s %10s %14s %10s\n','n','rhs','max rel err','verdict');
worst = 0;
for n = sizes
    A = randn(n); A = A*A' + n*eye(n);        % symmetric positive definite
    L = chol(A);                              % upper triangular, A = L'*L
    % gpml's infExact calls it with a vector and with eye(n); the MEX takes its
    % column count from the argument, so each shape is a separate code path.
    rhs = {'3 columns', randn(n,3); '1 column', randn(n,1); 'eye(n)', eye(n)};
    for k = 1:size(rhs,1)
        B = rhs{k,2};
        ref = A\B;
        try
            X = solve_chol(L,B);
        catch err
            fprintf(2,'  solve_chol threw: %s\n',err.message);
            return
        end
        relErr = max(abs(X(:)-ref(:)))/max(max(abs(ref(:))),1e-12);
        worst = max(worst,relErr);
        if relErr < 1e-10, verdict = 'ok'; else, verdict = '*** WRONG ***'; end
        fprintf('%6u %10s %14.3e %10s\n',n,rhs{k,1},relErr,verdict);
    end
end

fprintf('\nworst relative error: %.3e\n',worst);
if worst < 1e-10
    fprintf('The MEX agrees with MATLAB. solve_chol is NOT the problem.\n');
else
    fprintf(2,'\nThe MEX disagrees with MATLAB -- it is producing wrong answers.\n');
    fprintf(2,'This is enough to break the GP hyperparameter fits.\n');
    fprintf(2,'\nWorkaround: rename the MEX so the .m fallback is used instead:\n');
    d = fileparts(w);
    fprintf(2,'  movefile(''%s'',''%s'')\n', ...
        fullfile(d,['solve_chol.' mexext]),fullfile(d,['solve_chol.' mexext '.disabled']));
    fprintf(2,'The 148 gpml features then run on the MATLAB implementation --\n');
    fprintf(2,'same results, somewhat slower.\n');
end
fprintf('=======================================\n\n');
end
