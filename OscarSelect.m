
function [CoefMatrix dfMatrix SSMatrix] = OscarSelect(X, y, cvalues, propvalues, initcoef, method)

p = length(X(1,:));

% Standardize predictors to mean zero, variance 1, and center response to
% mean zero.

for i = 1:p
  X(:,i) = (X(:,i)-mean(X(:,i)))/std(X(:,i));
end;
y = y-mean(y);

if nargin < 6
    method = 2;
    if nargin < 5
        initcoef = [];
    end;
end;

cvalues=sort(cvalues);
propvalues=sort(propvalues);

if isempty(initcoef)
    [initcoef] = regress(y,X);     
end;
if length(initcoef)<p
    error('initial estimate must be of length p');
end;
if min(cvalues)<0
    error('all values of c must be nonnegative');
end; 
if max(cvalues)>1
    error('values for c cannot exceed 1');
end; 
if min(propvalues)<=0
    error('values for proportion must be greater than 0');
end;
if max(propvalues)>=1
    error('values for proportion must be smaller than 1');
end;

if method == 1
    [CoefMatrix dfMatrix SSMatrix] = OscarReg(X, y, cvalues, propvalues, initcoef);
end;
if method ~= 1
    [CoefMatrix dfMatrix SSMatrix] = OscarSeqReg(X, y, cvalues, propvalues, initcoef);
end;    