

% X should be a matrix whose rows are the observations and columns are the
% predictors (n by p). The intercept is omitted (so the response and each
% predictor should have mean zero).

function [CoefMatrix dfMatrix SSMatrix] = OscarSeqReg(X, y, cvalues, propvalues, initcoef)

p = length(X(1,:));
   
% Standardize predictors to mean zero, variance 1, and center response to
% mean zero.

for i = 1:p
  X(:,i) = (X(:,i)-mean(X(:,i)))/std(X(:,i));
end;
Xmatrix = [X -X]; % Need for positive and negative parts of beta
y = y-mean(y);

% Order initial estimate by magnitude, including ties.
% Initial estimate is used for first guess at ordering of coefficients and also to set maximum value
% for the bound t used in the constraint. The bound is given as proportion
% of the initial value.

[initcoeford, currorder] = sort(-abs(initcoef));
sameaslast = [0; (initcoeford(2:p) == initcoeford(1:(p-1)))];
startblocksame = [((sameaslast(2:p) - sameaslast(1:(p-1))) > 0); 0];
endblocksame = [((sameaslast(2:p) - sameaslast(1:(p-1))) < 0); sameaslast(p)];
nblocksame = sum(startblocksame);
vi = (1:p)';
visbs = vi(logical(startblocksame));
viebs = vi(logical(endblocksame));
for i = 1:nblocksame; 
    blockmean = mean(vi(visbs(i):viebs(i)));
    vi(visbs(i):viebs(i)) = blockmean * ones(viebs(i) - visbs(i) + 1,1);
end;
[tempinvsort,vind] = sort(currorder);
a1 = vi(vind)';
initcoeford = -initcoeford;   

% Grid search over c values, then by proportion of bound. For each c, the
% weighting is recomputed.

CoefMatrix = zeros(p,length(propvalues),length(cvalues));
SSMatrix = zeros(1,length(propvalues),length(cvalues));
dfMatrix = zeros(1,length(propvalues),length(cvalues));

for ccount = 1:length(cvalues)
    OrderMatrix = a1;
    weighting = a1;
    cvalue = cvalues(ccount);
    for i=1:p
        weighting(i) = (1-cvalue)+cvalue*(p-i);
    end;    
    for propcount = 1:length(propvalues)
        tbound = propvalues(propcount)*weighting*initcoeford;
        if (ccount == 1)
            if (propcount == 1)
                start=zeros(2*p,1);
            end;
        elseif (propcount > 1)
            start=[max(CoefMatrix(:,propcount-1,ccount),0);-min(CoefMatrix(:,propcount-1,ccount),0)];
        else
            start=[max(CoefMatrix(:,propcount,ccount-1),0);-min(CoefMatrix(:,propcount,ccount-1),0)];
        end;
        [coefs df ssquares conv] = OscarSeqOpt(tbound, cvalue, Xmatrix, y, start, OrderMatrix);
        CoefMatrix(:,propcount,ccount) = coefs;
        SSMatrix(:,propcount,ccount) = ssquares;
        dfMatrix(:,propcount,ccount) = df;
        if conv == 0
            fprintf('Optimization may not have converged properly for c = %g and prop = %g.\n', cvalue, propvalues(propcount));
        else
            fprintf('Optimization complete for c = %g and prop = %g.\n', cvalue, propvalues(propcount));
        end;
    end;
end;


