function groups = groupCandidateRegions(candidateRC, groupingRadiusPx)
    import membrane_tracking.focused_ism.*

    nCandidates = size(candidateRC, 1);
    groups = cell(0, 1);
    unused = true(nCandidates, 1);
    while any(unused)
        seed = find(unused, 1, 'first');
        queue = seed;
        unused(seed) = false;
        component = zeros(0, 1);

        while ~isempty(queue)
            current = queue(1);
            queue(1) = [];
            component(end+1, 1) = current; %#ok<AGROW>
            remaining = find(unused);
            if isempty(remaining)
                continue;
            end
            distance = hypot( ...
                candidateRC(remaining,1) - candidateRC(current,1), ...
                candidateRC(remaining,2) - candidateRC(current,2));
            neighbors = remaining(distance <= groupingRadiusPx);
            unused(neighbors) = false;
            queue = [queue; neighbors]; %#ok<AGROW>
        end
        groups{end+1, 1} = component; %#ok<AGROW>
    end
end
