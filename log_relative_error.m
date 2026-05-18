function LRE = log_relative_error(comp, exact)
    LRE = NaN(length(comp),1);
    for i = 1:length(comp)
        c = comp(i);
        [~,idx] = min(abs(exact - c));
        LRE(i) = -log10(abs(exact(idx) - c) / abs(exact(idx)));
        if isinf(LRE(i))
            LRE(i) = 16;
        end
    end
end