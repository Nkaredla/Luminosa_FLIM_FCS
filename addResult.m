function results = addResult(results, name, pass, detail)
%ADDRESULT Append one check to a test-result struct array.
%
% Kept as its own file rather than a local subfunction because this project
% requires one function per file (R2022a package prefixes make local helpers
% awkward to share between test scripts).

    results(end + 1) = struct('name', name, 'pass', logical(pass), ...
        'detail', detail);
end
