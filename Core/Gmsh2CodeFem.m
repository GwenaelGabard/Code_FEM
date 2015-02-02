function Gmsh2CodeFem(infile,outfile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads gmsh output file .msh and converts it to the CodeFEM format
% NODE ELEMENT ELEMENT_DOMAIN N_NODE N_ELEMENT...  
% Note: 
% Each domain is defined as multiples of 1000
% and each boundary by 1000*domain_i + boundary_i
% such that 4003 is the third boundary of domain number 4
%
% Elements supported:
% 1D: L2,T3
% 2D: Q4,T3,Q8,T6, T10, T15
% 3D: H8,T4,T10
%
% http://www.geuz.org/gmsh/doc/texinfo/gmsh_10.html#SEC62
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(infile,'rt');
if fid < 0
    error('file not found !!');
end
%
fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid);
%
N_NODE = str2double(fgetl(fid));
NUMNODE = zeros(1,N_NODE);
temp = fscanf(fid,'%f%f%f%f',[4,N_NODE]);
% -----------------
% COORDINATES
% -----------------
NODE(1:3,temp(1,:)) = [temp(2,:); temp(3,:); temp(4,:)];

% ------------------
% GMSH NUMEROTATION
% ------------------
NUMNODE(temp(1,:))=1:N_NODE; % Gmesh numerotation index
%
%%%%%%%%%%%%%%%%%%%%%  ELEMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = fgetl(fid); %
temp = fgetl(fid); % on lit la ligne $ENDNOD
temp = fgetl(fid); % on lit la ligne $ELEM

% 
N_ELEMENT = str2double(fgetl(fid));
ELEMENT = zeros(1,N_ELEMENT);
ELEMENT_DOMAIN = zeros(3,N_ELEMENT);
N_NODE_ELEMENT = zeros(1,N_ELEMENT);
NODE_DOMAIN = ones(1,N_NODE);
N_DIM = 2;

%
for k=1:N_ELEMENT
    %
    out = str2num(fgetl(fid));
    type = out(2);
    
    if out(3) ~= 2
        error('Physical groups should be defined !!!')
    end
    % pointeur debut liste des noeuds
    point = 6;
    
    %------------------------
    % WHICH DOMAIN/BOUNDARY
    %------------------------
    %
    R = mod(out(4),1000);
    Q = (out(4) - R) / 1000;
    if R == 0
        ELEMENT_DOMAIN(1:3,k) = [Q ; 0 ; 1];
    else
        ELEMENT_DOMAIN(1:3,k) = [Q ; 1; R];
    end
    ELEMENT_DOMAIN(4,k) = type;
    
    %-------------------
    % CONNECTIVITY
    %-------------------
    %
    if type==1
        %**************************************
        %           LINEAR LINE
        %**************************************`
        ELEMENT(1:2,k) = NUMNODE(out(point:end)).';
        N_NODE_ELEMENT(k) = 2;
    elseif type==3
        %**************************************
        %            LINEAR QUADRANGLE
        %**************************************`
        ELEMENT(1:4,k) = NUMNODE(out(point:end)).';
        N_NODE_ELEMENT(k) = 4;
    elseif type==2
        %**************************************
        %            LINEAR TRIANGLE
        %**************************************
        ELEMENT(1:3,k) = NUMNODE(out(point:end)).';
        N_NODE_ELEMENT(k) = 3;
    elseif type==8
        %**************************************
        %           QUADRATIC LINE
        %**************************************
        temp = out(point:end);
        ELEMENT(1:3,k) = NUMNODE(temp([1 3 2])).';
        N_NODE_ELEMENT(k) = 3;
        ELEMENT_DOMAIN(1:3,k) = [ Q 1 out(4)-(Q*1000) ];
    elseif type==10
        %**************************************
        %           QUADRATIC QUADRANGLE
        %**************************************
        temp = out(point:end);
        ELEMENT(1:9,k) = NUMNODE(temp([1 5 2 6 3 7 4 8 9])).';
        N_NODE_ELEMENT(k) = 9;
        ELEMENT_DOMAIN(1:3,k) = [ Q 0 1 ];
    elseif type==9
        %**************************************
        %           QUADRATIC TRIANGLE
        %**************************************
        temp = out(point:end);
        ELEMENT(1:6,k) = NUMNODE(temp([1 4 2 5 3 6])).';
        N_NODE_ELEMENT(k) = 6;
        ELEMENT_DOMAIN(1:3,k) = [ Q 0 1 ];
    elseif type==4
        %**************************************
        %           LINEAR TETRAHEDRON
        %**************************************
        ELEMENT(1:4,k) = NUMNODE(out(point:end)).';
        N_NODE_ELEMENT(k) = 4;
        N_DIM = 3;
    elseif type==5
        %**************************************
        %           LINEAR HEXAHEDRON
        %**************************************
        ELEMENT(1:8,k) = NUMNODE(out(point:end)).';
        N_NODE_ELEMENT(k) = 8;
        N_DIM = 3;
    elseif type==11
        %**************************************
        %           QUADRATIC TETRAHEDRON
        %**************************************
        temp = out(point:end);
        ELEMENT(1:10,k) = NUMNODE(temp([1 5 2 6 3 7 8 10 9 4])).';
        N_NODE_ELEMENT(k) = 10;
        N_DIM = 3;
    elseif type==26
        %**************************************
        %           CUBIC LINE
        %**************************************
        temp = out(point:end);
        ELEMENT(1:4,k) = NUMNODE(temp([1 3 4 2])).';
        N_NODE_ELEMENT(k) = 4;
        %ELEMENT_DOMAIN(1:3,k) = [ Q 1 out(4)-(Q*1000) ]; ?????
    elseif type==21
        %**************************************
        %           CUBIC TRIANGLE
        %**************************************
        temp = out(point:end);
        ELEMENT(1:10,k) = NUMNODE(temp([1 4 5 2 6 7 3 8 9 10])).';
        N_NODE_ELEMENT(k) = 10;
        %ELEMENT_DOMAIN(1:3,k) = [ Q 0 1 ];????
   elseif type==36
        %**************************************
        %           CUBIC QUADRANGLE
        %**************************************
        temp = out(point:end);
%         ELEMENT(1:16,k) = NUMNODE(temp([1 5 6 2 7 8 3 9 10 4 11 12 14 15 16 13])).';
        ELEMENT(1:16,k) = NUMNODE(temp([1 5 6 2 7 8 3 9 10 4 11 12 13 14 15 16])).';
        N_NODE_ELEMENT(k) = 16;
        %ELEMENT_DOMAIN(1:3,k) = [ Q 0 1 ];????
     elseif type==27
        %**************************************
        %           QUARTIC LINE
        %**************************************
        temp = out(point:end);
        ELEMENT(1:5,k) = NUMNODE(temp([1 3 4 5 2])).';
        N_NODE_ELEMENT(k) = 5;
        %ELEMENT_DOMAIN(1:3,k) = [ Q 1 out(4)-(Q*1000) ]; ?????
    elseif type==23
        %**************************************
        %           QUARTIC TRIANGLE
        %**************************************
        temp = out(point:end);
        ELEMENT(1:15,k) = NUMNODE(temp([1 4 5 6 2 7 8 9 3 10 11 12 13 14 15])).';
        N_NODE_ELEMENT(k) = 15;
        %ELEMENT_DOMAIN(1:3,k) = [ Q 0 1 ];????
    elseif type==37
        %**************************************
        %           QUARTIC QUADRANGLE
        %**************************************
        temp = out(point:end);
%         ELEMENT(1:25,k) = NUMNODE(temp([1 5 6 7 2 8 9 10 3 11 12 13 4 14 15 16 18 22 19 23 20 24 17 21 25])).';
        ELEMENT(1:25,k) = NUMNODE(temp([1 5 6 7 2 8 9 10 3 11 12 13 4 14 15 16 17 21 18 22 19 23 20 24 25])).';
        N_NODE_ELEMENT(k) = 25;
        %ELEMENT_DOMAIN(1:3,k) = [ Q 0 1 ];????
    end
end

%
fclose(fid);

ELEMENT_ID = 1:N_ELEMENT;
NODE_ID = 1:N_NODE;

% détermination de N_DOMAIN
N_DOMAIN = max(ELEMENT_DOMAIN(1,:));

for k = N_DOMAIN:-1:1
    Elt_list = find(ELEMENT_DOMAIN(1,:) == k);
    Node_list = ELEMENT(:,Elt_list);
    Node_list = Node_list(:);
    Node_list = Node_list(find(Node_list ~= 0));
    NODE_DOMAIN(1,Node_list') = k*ones(size(Node_list'));
end

%if N_DIM == 2
%    NODE = NODE(1:2,:);
%end


% Sauvegarde du maillage dans un fichier Code FEM
save(outfile,'N_DIM','N_NODE','NODE','N_ELEMENT','N_NODE_ELEMENT','ELEMENT','N_DOMAIN','ELEMENT_DOMAIN','NODE_DOMAIN','ELEMENT_ID','NODE_ID');
