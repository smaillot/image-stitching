function H = solveHomo(match1,match2)

    n = length(match1);
    A = zeros(2*n,9);
    A(:,1:3) = kron([match1 ones(n,1)], [1;0]);
    A(:,4:6) = kron([match1 ones(n,1)], [0;1]);
    x1 = match1(:,1)';
    y1 = match1(:,2)';
    x2 = match2(:,1)';
    y2 = match2(:,2)';
    A(1:2:2*n,7) = -x2.*x1;
    A(2:2:2*n,7) = -y2.*x1;
    A(1:2:2*n,8) = -x2.*y1;
    A(2:2:2*n,8) = -y2.*y1;
    A(1:2:2*n,9) = -x2;
    A(2:2:2*n,9) = -y2;

    [A,B] = eig(A'*A);
    H = reshape(A(:,1),[3,3])';
    H = H/H(end);

end