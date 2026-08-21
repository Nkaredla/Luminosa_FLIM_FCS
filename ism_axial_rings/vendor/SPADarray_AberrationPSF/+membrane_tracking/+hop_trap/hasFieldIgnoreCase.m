function tf = hasFieldIgnoreCase(s, name)
%HASFIELDIGNORECASE Test a structure field name without case sensitivity.

    tf = any(strcmpi(fieldnames(s), name));
end
