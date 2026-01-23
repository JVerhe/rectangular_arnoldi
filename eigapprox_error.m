function err = eigapprox_error(exact,comp,shift)

    comp = comp(~isnan(comp));
    comp = comp(~isinf(comp));

    err = 0;
    for i = 1:length(comp)
        if isempty(exact)
            err = err + sum(abs(comp.^2));  
            break
        else
            [match,idx] = min(exact - comp(i));
            err = err + abs(match^2);
            % comp(1) = [];
            % exact(idx) = [];
        end
    end
    err = sqrt(err);
end