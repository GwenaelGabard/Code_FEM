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


% Provide a simple GUI to select what is shown on screen

% The prefix of the post-processing functions
FILE_ID = 'POST_';

namelist = cell(N_DOMAIN,1);
domainname = cell(N_DOMAIN+1,1);
filename = cell(N_DOMAIN+1,1);
for k=1:N_DOMAIN
   domainname{k} = ['Domain' num2str(k)];
   addpath(DOMAIN_PATH{k})
   names = dir(DOMAIN_PATH{k});
   n = size(names);
   number = 0;
   for j=1:n
      if (~names(j).isdir)&strncmp(names(j).name,FILE_ID,size(FILE_ID,2))
         number = number+1;
         LIST(number) = j;
      end
   end
   namelist{k} = cell(number+2,1);
   filename{k} = cell(number,1);
   for j=1:number
      temp = names(LIST(j)).name;
      filename{k}{j} = temp(1:end-2);
      namelist{k}{j} = feval(temp(1:end-2));
   end
   namelist{k}{number+1} = 'BACK';
   namelist{k}{number+2} = 'END';
   rmpath(DOMAIN_PATH{k})
end
domainname{N_DOMAIN+1} = 'END';

A = 1;
while A==1
   a = menu('Post Processing',domainname);
   if a==N_DOMAIN+1
      A = 0;
   else
      addpath(DOMAIN_PATH{a})
      B = 1;
      while B==1
         b = menu(domainname{a},namelist{a});
         if b==length(namelist{a})
            A = 0;
            B = 0;
         end
         if b==length(namelist{a})-1
            B = 0;
         end
         if b<length(namelist{a})-1
            feval(filename{a}{b},a);
         end
      end
      rmpath(DOMAIN_PATH{a})
   end
end

% Tidy up the workspace
clear k j number temp A B a b domainname filename namelist n names
