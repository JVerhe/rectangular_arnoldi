function [B_new, P] = align_and_minimize_tol(A, B, tol)
    if nargin < 3
        tol = 1e-10; 
    end
    B_new = zeros(size(A));
    P = zeros(size(A));
    
    matched_A_idx = false(size(A));
    matched_B_idx = false(size(B));
    
    for i = 1:length(A)
        idx = find(abs(B - A(i)) <= tol & ~matched_B_idx, 1);
        
        if ~isempty(idx)
            B_new(i) = B(idx);
            P(i) = idx;
            matched_A_idx(i) = true;
            matched_B_idx(idx) = true;
        end
    end
    
    open_positions = find(~matched_A_idx);
    A_open = A(open_positions);
    
    B_rem_orig_indices = find(~matched_B_idx);
    B_rem = B(B_rem_orig_indices);
    
    [~, sort_idx_A] = sort(A_open);
    [B_rem_sorted, sort_idx_B] = sort(B_rem);
    
    B_new(open_positions(sort_idx_A)) = B_rem_sorted;
    
    P(open_positions(sort_idx_A)) = B_rem_orig_indices(sort_idx_B);
end