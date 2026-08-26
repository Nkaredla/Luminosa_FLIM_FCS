function mark = ternaryPassMark(pass)
%TERNARYPASSMARK 'PASS' or 'FAIL' for one logical, for aligned test output.

    if pass
        mark = 'PASS';
    else
        mark = 'FAIL';
    end
end
