function results = run_all_tests()
%RUN_ALL_TESTS Execute the dynamic APV-AD unit test suite.
%
%   RESULTS = RUN_ALL_TESTS() runs APVADTest and prints a summary.
%
%   This is the MATLAB counterpart of the project's pytest suite. Run it
%   after changing any sub-model, and in particular after changing a
%   parameter that several models share -- the tests are written against
%   closed forms and identities rather than stored numbers, so they remain
%   valid when parameters move and will still fail if the physics breaks.
%
%   See also APVADTEST, SMOKE_TEST.

here = fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(here)));

results = runtests('APVADTest');

fprintf('\n================ test summary ================\n');
fprintf('  passed   : %d\n', nnz([results.Passed]));
fprintf('  failed   : %d\n', nnz([results.Failed]));
fprintf('  incomplete: %d\n', nnz([results.Incomplete]));
fprintf('  duration : %.1f s\n', sum([results.Duration]));

failed = results([results.Failed]);
for k = 1:numel(failed)
    fprintf('  FAILED: %s\n', failed(k).Name);
end
fprintf('==============================================\n');
end
