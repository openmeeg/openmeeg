function [data] = load_bin_sym(filename)

% LOAD_TXT_SYM   Load binary symmetric Matrix
%
%   Load binary symmetric Matrix
%   
%   SYNTAX
%       [DATA] = LOAD_BIN_SYM(FILENAME)
%   

%
%   Created by Alexandre Gramfort on 2007-09-26.
%

file = fopen(filename,'r');
dim = fread(file,1,'uint32');
data = zeros(dim,dim);
data(tril(ones(dim,dim)) > 0) = fread(file,dim*(dim+1)/2,'double');
data = data + data' - diag(diag(data));
fclose(file);

