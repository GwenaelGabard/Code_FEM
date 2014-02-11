% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Copyright (c) 2014 Code_FEM developers


function [number,list]=ADD_TO_LIST(list,token)
% list : a list of strings
% token : a string
% If token is already in the list then list remain unchanged
% and number is the index where token is in the list.
% If token isn't in the list, we add it to the list and return the
% corresponding index in the list in number.
% Upper and lower cases are different.
% Example :
% list = {'alpha','beta','Alpha'};
% [number,list] = ADD_TO_LIST(list,'beta')
% [number,list] = ADD_TO_LIST(list,'Beta')

temp = find(strcmp(list,token));
if isempty(temp)
    % token is not yet in the list
    list = [list cellstr(token)];
    number = length(list);
else
    % token is already in the list
    number = temp;
end
