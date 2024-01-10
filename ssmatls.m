function [y, Qmat, Rmat, h]=ssmatls(n);
% ssmatls smoothing matrix for smoothing spline method
%    the definition is based on page 12-13 in Silverman (1994)
%
% (c) copyright Lingsong Zhang (lingsong@purdue.edu)

h=1/(n+1);

%
Qmat=sparse(n, n-1);
for (j=2:n-1);
    Qmat(j-1, j)=1/h;
    Qmat(j, j)=-2/h;
    Qmat(j+1, j)=1/h;
end;
Qmat=Qmat(:, 2:(n-1));

Rmat=sparse(n-1, n-1);
for (i=1:(n-2));
    Rmat(i, i)=(2/3)*h;
    Rmat(i, i+1)=h/6;
    Rmat(i+1, i)=h/6;
end;
Rmat(n-1, n-1)=(2/3)*h;
Rmat=Rmat(2:(n-1), 2:(n-1));

y=Qmat*inv(Rmat)*Qmat';
