function [data] = load_txt_sym(filename)

% LOAD_TXT_SYM   Load ascii symmetric matrix
%
%   Load ascii symmetric matrix
%
%   SYNTAX
%       [DATA] = LOAD_TXT_SYM(FILENAME)
%

%
%   Created by Alexandre Gramfort on 2007-09-26.
%

file = fopen(filename);
rawdata = textscan(file,'%f');
rawdata = cell2mat(rawdata);
dim = (-1 + sqrt(1+8*length(rawdata)))/2;
assert(dim == ceil(dim),'Bad dimension for a symmetric matrix')
data = zeros(dim,dim);
data(tril(ones(dim,dim)) > 0) = rawdata;
data = data + data' - diag(diag(data));
fclose(file);

