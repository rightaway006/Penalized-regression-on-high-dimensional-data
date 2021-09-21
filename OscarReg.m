

% X should be a matrix whose rows are the observations and columns are the
% predictors (n by p). The intercept is omitted (so the response and each
% predictor should have mean zero).

function [CoefMatrix dfMatrix SSMatrix] = OscarReg(X, y, cvalues, propvalues, initcoef)

p = length(X(1,:));
q=p*(p-1)/2;


% Standardize predictors to mean zero, variance 1, and center response to
% mean zero.

for i = 1:p
  X(:,i) = (X(:,i)-mean(X(:,i)))/std(X(:,i));
end;
y = y-mean(y);



corrmat=X'*X;
[ivec,jvec,svec]=find(corrmat);
Sparseavec=[ivec;ivec+p;ivec;ivec+p];
Sparsebvec=[jvec;jvec+p;jvec+p;jvec];
Element1vec=[svec;svec;-svec;-svec];

F=sparse(Sparseavec,Sparsebvec,Element1vec,2*p+q,2*p+q);

clear corrmat ivec jvec svec Sparseavec Sparsebvec Element1vec;

c=[-X'*y;X'*y;zeros(q,1)];
lowbound=[zeros(2*p,1);-inf(q,1)];

Sparse1vec=[ones(2*p+q,1);2;2;2;3;3;3];
Sparse2vec=[(1:(2*p+q))';1;p+1;2*p+1;2;p+2;2*p+1];
Elementvec=[ones(2*p,1);ones(q,1);1;1;-1;1;1;-1];
rowcount=2;
paircount=2*p+1;

initcoef1=[max(initcoef,0);-min(initcoef,0)];

for i=1:p
    for j=(i+1):p
        initcoef1=[initcoef1;max(abs(initcoef(i)),abs(initcoef(j)))];
        if rowcount>2
            Sparse1vec=[Sparse1vec;rowcount;rowcount;rowcount;rowcount+1;rowcount+1;rowcount+1];
            Sparse2vec=[Sparse2vec;i;i+p;paircount;j;j+p;paircount];
            Elementvec=[Elementvec;1;1;-1;1;1;-1];
        end;
        rowcount=rowcount+2;
        paircount=paircount+1;
    end;
end;


% Grid search over c values, then by proportion of bound. For each c, the
% constraint matrix is updated.

CoefMatrix = zeros(p,length(propvalues),length(cvalues));
SSMatrix = zeros(1,length(propvalues),length(cvalues));
dfMatrix = zeros(1,length(propvalues),length(cvalues));

for ccount = 1:length(cvalues)
    cvalue = cvalues(ccount);
    tempvec=[(1-cvalue)*ones(2*p,1);cvalue*ones(q,1)];
    Elementvec(1:(2*p+q))=tempvec;
    clear tempvec;

    A=sparse(Sparse1vec,Sparse2vec,Elementvec,2*q+1,2*p+q);
    initialnorm=A(1,:)*initcoef1; 
    for propcount = 1:length(propvalues)
        tbound = propvalues(propcount)*initialnorm;
        b=[tbound; zeros(2*q,1)];
        
        if (ccount == 1)
            if (propcount == 1)
                start=zeros(2*p+q,1);
            end;
        elseif (propcount > 1)
            start=[max(CoefMatrix(:,propcount-1,ccount),0);-min(CoefMatrix(:,propcount-1,ccount),0)];
            for i=1:p
                for j=(i+1):p
                    start=[start;max(abs(start(i)),abs(start(j)))];
                end;
            end;
        else
            start=[max(CoefMatrix(:,propcount,ccount-1),0);-min(CoefMatrix(:,propcount,ccount-1),0)];
            for i=1:p
                for j=(i+1):p
                    start=[start;max(abs(start(i)),abs(start(j)))];
                end;
            end;
        end;

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%
%%%%%%
%%%%%% The TOMLAB package is now used as a quadratic programming solver. 

%%%%%% IF AN ALTERNATIVE QUADRATIC PROGRAMMING SOLVER IS AVAILABLE, 
%%%%%% REPLACE THIS SECTION WITH A CALL TO THAT SOLVER BASED ON  
%%%%%% THE FORMULATION 

% minimize (over x) : 0.5 x' F x + c' x
% subject to A x <= b and x >= lowbound
%  

        Prob = qpAssign(F, c, A, [], b, lowbound, [], start, [], [], [], [], [], [], [], []);
        Result = tomRun('sqopt', Prob, 0);

        x=Result.x_k;
        exitflag = Result.ExitFlag;
        conv = (exitflag == 0);
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
        CoefMatrix(:,propcount,ccount) = round((x(1:p)-x(p+1:2*p))*10^7)*10^(-7);
        SSMatrix(:,propcount,ccount) = sumsqr(y-X*CoefMatrix(:,propcount,ccount));
    
        EffParVec=unique(abs(CoefMatrix(:,propcount,ccount)));
        EffParVec=EffParVec(EffParVec>0);
        dfMatrix(:,propcount,ccount) = length(EffParVec);
    
        if conv == 0
            fprintf('Optimization may not have converged properly for c = %g and prop = %g.\n', cvalue, propvalues(propcount));
        else
            fprintf('Optimization complete for c = %g and prop = %g.\n', cvalue, propvalues(propcount));
        end;
        
    end;
    
end;


