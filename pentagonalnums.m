function PTK = pentagonalnums(m)
% Computes the first m generalized pentagonal numbers and returns an array of them
% See https://en.wikipedia.org/wiki/Pentagonal_number
% Used to calculate the partition function of an integer as in 
% https://en.wikipedia.org/wiki/Partition_(number_theory)
% Called by partitionnum.m
% Jason D. Kahn 11/30/18
%
PTK = zeros(m,1);
    if mod(m,2) == 0
        for ll = 1:m/2
            PTK(2*ll-1) = ll*(3*ll-1)/2;
            PTK(2*ll) = -ll*(3*-ll-1)/2;
        end
    else
        for ll = 1:(m-1)/2
            PTK(2*ll-1) = ll*(3*ll-1)/2;
            PTK(2*ll) = -ll*(3*-ll-1)/2;
        end
        ll = (m-1)/2+1;
        PTK(2*ll-1) = ll*(3*ll-1)/2;
    end
end