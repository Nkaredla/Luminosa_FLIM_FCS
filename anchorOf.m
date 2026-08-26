function value = anchorOf(ref)
%ANCHOROF The SLB lifetime a reference measurement settled on.
%
% The pooled bare-SLB fit may be a mono fit (field tauNs) or the two-component
% fit whose SHORT lifetime is the anchor (field slbTauNs). This picks whichever
% is present so callers do not have to know which model was used.

    if isfield(ref, 'anchorNs') && ~isempty(ref.anchorNs)
        value = ref.anchorNs;
    elseif isfield(ref.pooled, 'slbTauNs')
        value = ref.pooled.slbTauNs;
    else
        value = ref.pooled.tauNs;
    end
end
