function index = binaryBracketIncreasing( A, num )
%--------------------------------------------------------------------------
% Syntax:       [index] = binarySearch(A, n, num);
%               
% Inputs:       A: Array (sorted, asc) that you want to search
%               num: Number you want to search in array A
%               
% Output:       index: Return position in A that A(index) <= num < A(index+1)
%                      or -1 if num < A(1)
%               if num > A(end), return size(A)
%               
% Description:  Use binary search to find left bracket for increasing function
%               
% Complexity:   O(1)    best-case performance
%               O(log_2 (n))    worst-case performance
%               O(1)      auxiliary space
%--------------------------------------------------------------------------
n = length(A);
left = 1;
right = n;

if num <= A
    index = 1;
elseif num >= A
    index = n;
else
    for i = 1:n
        mid = ceil((left + right) / 2);

        if A(mid) <= num
            if A(mid+1) > num
                index = mid;
                break;
            else
                left = mid;
            end
        else
            right = mid;
        end
    end

end

