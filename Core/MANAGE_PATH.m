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

global LIBRARY_PATH

DIRECTORY = fullfile('Code_FEM','Core');

% Get the directories in the path
P = path;
numbers = find(P==pathsep);
for k=0:length(numbers)
    if k==0
        token = P(1:(numbers(1)-1));
    end
    if k==length(numbers)
        token = P((numbers(k)+1):end);
    end
    if (k>0)&&(k<length(numbers))
        token = P((numbers(k)+1):(numbers(k+1)-1));
    end
    if strcmpi(token(end-length(DIRECTORY)+1:end),DIRECTORY)
        break
    end
end

if ~strcmpi(token(end-length(DIRECTORY)+1:end),DIRECTORY)
    error('Can''t find the directory');
else
    LIBRARY_PATH{1} = token(1:end-5);
end



DIRECTORY_PRIVATE = fullfile('Code_FEM_Private','Core');
% Get the directories in the path
P = path;
numbers = find(P==pathsep);
for k=0:length(numbers)
    if k==0
        token = P(1:(numbers(1)-1));
    end
    if k==length(numbers)
        token = P((numbers(k)+1):end);
    end
    if (k>0)&&(k<length(numbers))
        token = P((numbers(k)+1):(numbers(k+1)-1));
    end
    if strcmpi(token(end-length(DIRECTORY_PRIVATE)+1:end),DIRECTORY_PRIVATE)
        break
    end
end

if strcmpi(token(end-length(DIRECTORY_PRIVATE)+1:end),DIRECTORY_PRIVATE)
    LIBRARY_PATH{2} = token(1:end-5);
end

DIRECTORY_LMS_Private = fullfile('Code_FEM_LMS_Private','Core');
% Get the directories in the path
P = path;
numbers = find(P==pathsep);
for k=0:length(numbers)
    if k==0
        token = P(1:(numbers(1)-1));
    end
    if k==length(numbers)
        token = P((numbers(k)+1):end);
    end
    if (k>0)&&(k<length(numbers))
        token = P((numbers(k)+1):(numbers(k+1)-1));
    end
    if strcmpi(token(end-length(DIRECTORY_LMS_Private)+1:end),DIRECTORY_LMS_Private)
        break
    end
end

if strcmpi(token(end-length(DIRECTORY_LMS_Private)+1:end),DIRECTORY_LMS_Private)
    LIBRARY_PATH{3} = token(1:end-5);
end

% Remove temporary variables
clear P k numbers token
