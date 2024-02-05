function QNK = partitionfct(n,k)
%Calculate the partition function, the number of ways to express an
% integer n as a sum of k terms. Uses a recurrence relation from
%   http://mathworld.wolfram.com/PartitionFunctionq.html
% (as of 2024 this page seems to be unreadable)
%  This file includes partitionfctrec, a recursive version of itself,
%   as a nested function so it has access to the workspace of the main
%   function. 
%  Calls partitionnum.m 
if (n < 0 || k<0)
    disp('Partition function undefined for negative args')
    return
end
Qmat(1,1) = 1;
Pmat(1) = 1;
if (k == 0)
	QNK = 0;
elseif (n == 1)
    QNK = 1;
elseif (k >= n)
    if (numel(Pmat) > n && Pmat(n) ~= 0)
        QNK = Pmat(n);
    else
        QNK = partitionnum(n);    
        Pmat(n) = QNK;
        Qmat(n,k) = QNK;
    end
else
    if (size(Qmat,1) > n && size(Qmat,2) > k && Qmat(n,k) ~= 0)
        QNK = Qmat(n,k);
    else
    %n
    %k
        QNK = partitionfctrec(n,k-1) + partitionfctrec(n-k,k);
        Qmat(n,k) = QNK;
    end
end
    
function QNK = partitionfctrec(m,l)
% This is the recursive function, should rarely be called directly
% It is defined here as a nested function so it has access to the calling
% program's workspace
    if (m < 0 || l < 0)
        disp('Partition function undefined for negative args')
        return
    end
    if (l == 0)
        QNK = 0;
    elseif (m == 1)
        QNK = 1;
    elseif (l >= m)
        if (numel(Pmat) > m) && (Pmat(m) ~= 0)
            QNK = Pmat(m);
        else
            QNK = partitionnum(m);    
            Pmat(m) = QNK;
            Qmat(m,l) = QNK;
        end
    else
        if (size(Qmat,1) > m && size(Qmat,2) > l && Qmat(m,l) ~= 0)
            QNK = Qmat(m,l);
        else
            QNK = partitionfctrec(m,l-1) + partitionfctrec(m-l,l);
            Qmat(m,l) = QNK;
        end
    end
end

end


