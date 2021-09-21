% Read Soil Dataset
data = xlsread('C:\Users\sarmad\Desktop\soil.xlsx');

% -------------- Separating features in separate variable ---------------
X = data(:, 1:15);
  
% -------------- Separating magnitude in separate variable --------------
y = data(:,16);




% Choose grid of parameter values to use
cvalues = [0; .01; .05; .1;.25;.5;.75;.9;1];
propvalues = [.0001; .001; .002; .00225; .0025; .00275; .0028; .0029; .003; .004; .005; .0075; .01;.025; .05; .1; .15;.2;.3;.4;.5;.6];

%%%% method = 2, chooses the sequential algorithm.

method = 2;     
initcoef = [];
[CoefMatrix dfMatrix SSMatrix] = OscarSelect(X, y, cvalues, propvalues, initcoef, method); % Calls function to perform optimization on the grid.

