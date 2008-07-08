function [data] = load_bin(filename)

% LOAD_BIN   Load binary Matrix
%
%   Load binary Matrix
%
%   SYNTAX
%       [DATA] = LOAD_BIN(FILENAME)
%

%
%   Created by Alexandre Gramfort on 2007-05-22.
%

file = fopen(filename,'r');
dims = fread(file,2,'uint32');
data = fread(file,prod(dims),'double');
data = reshape(data,dims');
fclose(file);
