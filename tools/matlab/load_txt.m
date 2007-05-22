function [data] = load_txt(filename)

% LOAD_TXT   Load ascii matrix
%
%   Load ascii matrix
%
%   SYNTAX
%       [DATA] = LOAD_TXT(FILENAME)
%

%
%   Created by Alexandre Gramfort on 2007-05-22.
%

data = load(filename);
