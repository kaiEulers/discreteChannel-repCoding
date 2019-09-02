% p2 is the probability that a bit error will occur
% (bit swapped)

function [x, y] = disChn_rep(p2)
    
    x = 0;
    yProb = rand;
    
    if (yProb <= (1 - p2))
        % (1 - p2) chance that output will be the same as input
        y = 0;
    else
        % p2 chance that output will be not be the same as
        % input
        y = 1;
    end
    
end