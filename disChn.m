% p2 is the probability that a bit error will occur
% (bit swapped)

function [x, y] = disChn(p2)
    
    xProb = rand;
    
    if (xProb < 0.5)
        
        x = 0;
        yProb = rand;
        
        if (yProb <= (1 - p2))
            y = 0;
        else
            y = 1;
        end
        
    else
        
        x = 1;
        yProb = rand;
        
        if (yProb <= (1 - p2))
            y = 1;
        else
            y = 0;
        end
        
    end
end