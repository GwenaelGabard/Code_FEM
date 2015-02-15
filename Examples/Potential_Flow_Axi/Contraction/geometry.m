function [x,y] = geometry(bs,s,x)
% CONVENTION AXI : R=X
%                  Z=Y

nbs = 6;
H = 1;
L = 1;
D = 2;
C = 0.3;

% Return the number of boundary segments
if nargin == 0,
   x = nbs;   
   return
end

d = [0 0 0 0 0 0 0 0 % start parameter value
     1 1 1 1 1 1 1 1 % end parameter value
     1 1 1 1 1 1 1 1 % left hand region
     0 0 0 0 0 0 0 0 % right hand region
    ];
bs1 = bs(:)';
if find(bs1<1 | bs1>nbs),
   error('Non existent boundary segment number')
end

if nargin == 1,
   x = d(:,bs1);
   return
end

if nargin == 2
   x = zeros(size(s));
   y = zeros(size(s));
   [m,n] = size(bs);
   if m==1 && n==1,
      bs = bs*ones(size(s)); % expand bs
   elseif m~=size(s,1) || n~=size(s,2),
      error('bs must be scalar or of same size as s');
   end

   if ~isempty(s)
      ii = find(bs==1);
      x(ii) = 0 + (L+D+L)*s(ii);
      y(ii) = 0;
      ii = find(bs==2);
      x(ii) = L+D+L;
      y(ii) = H*s(ii);
      ii = find(bs==3);
      x(ii) = L+D+L - L*s(ii);
      y(ii) = H;
      ii = find(bs==4);
      x(ii) = L+D - D*s(ii);
      y(ii) = H - C*(1-cos(2*pi*s(ii)))/2;
      ii = find(bs==5);
      x(ii) = L - L*s(ii);
      y(ii) = H;
      ii = find(bs==6);
      x(ii) = 0;
      y(ii) = H*(1-s(ii));
   end
   return
end
