function shuf = segperm(M, dim, n)
    
    % Preprocessing
    if dim == 1
        M = M';
    elseif dim == 3
        fprintf("Nope! Try a different sized array");
    else 
        M = M;
    end

    % If n is a proportion
    if n < 1
        n = n * b;
    end
    
    % Get size
    [a, b] = size(M);

    % Segment
    segs = 1:round(n):b;
    b = 1:b;
    
    % Randperm
    idx = zeros(size(M));
    for i = 1:a
        id = [];
        p = randperm(length(segs));
        for j = 1:length(segs)
            low = segs(p(j));
            upp = segs(p(j)) + round(n);
            idb = find(b >= low & b < upp);
            id = cat(2, id, idb);
        end
        idx(i,:) = id;
    end

    % Shuffled Array
    shuf = M(idx);
end