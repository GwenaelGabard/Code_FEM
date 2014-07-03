function [Ke,Fe,Re,Ve] = DBEM_Wall_Wall(m,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helmholtz equation BEM formulation (direct BEM)
% Isoparametric high-order element (Lagrange polynomial)
% Interpolation depends on order
% L2->1, L3->2, ...
% Right hand side: Diffraction by a plane wave
% Formulation:
% 1/2*\int_Sx {\phi*(X) \phi(Y)} dSx
% - \int_Sx \int_Sy {\phi*(X) \phi(Y) dG(X,Y)/dny} dSy dSx
% = \int_Sx {\phi*(X) \phi_inc(X)} dSp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Order

NrOfNodes = (Order+1);

if nargin==0
    Ke = [1 NrOfNodes ones(1,NrOfNodes) 0];
    return
end

global omega 
global ELEMENT NODE NODE_DATA
global N_GAUSS_POINT GAUSS_POINT GAUSS_WEIGHT 
global SHAPEFUNC DSHAPEFUNC

RHO = NODE_DATA(1,ELEMENT(1,m));
C0 = NODE_DATA(2,ELEMENT(1,m));
K = (omega/C0);

Ke = zeros(NrOfNodes,NrOfNodes);
Fe = zeros(NrOfNodes,1);
Re = zeros(NrOfNodes,NrOfNodes);
Ve = zeros(NrOfNodes,1);

% Density quadrature factor (For the auto-influence computation)
QuadAutoFactor = 1;

x_nodem = NODE(1,ELEMENT(:,m))';
y_nodem = NODE(2,ELEMENT(:,m))';

%--------------------------------------------------------------------------
%     Double integration (Calling function with a double argument)
%--------------------------------------------------------------------------
if nargin ==2
    x_nodep = NODE(1,ELEMENT(:,p))';
    y_nodep = NODE(2,ELEMENT(:,p))';
    % Standard quadrature (elements do not coincide)
    if m~=p
        % integrating over M
        for mm=1:N_GAUSS_POINT
            u = GAUSS_POINT(mm);
            Nu = SHAPEFUNC(:,mm)';
            dNu = DSHAPEFUNC(:,mm)';
            Jm = [dNu*x_nodem dNu*y_nodem];
            M = [Nu*x_nodem Nu*y_nodem];
            % integrating over P
            for pp=1:N_GAUSS_POINT
                v = GAUSS_POINT(pp);
                Nv = SHAPEFUNC(:,pp)';
                dNv = DSHAPEFUNC(:,pp)';
                Jp = [dNv*x_nodep dNv*y_nodep];
                P = [Nv*x_nodep Nv*y_nodep];
                R = norm(P-M);
                %G = -i/4*besselh(0,1,K*R);
                np = [-Jp(2) ; Jp(1)]/norm([-Jp(2) ; Jp(1)]);
                dGdnp = (1i*K/4)*besselh(1,1,K*R)*(P-M)/R*np;
                J = norm(Jm)*norm(Jp);
                % Compute the operator
                Ke = Ke - GAUSS_WEIGHT(mm)*GAUSS_WEIGHT(pp)*Nu'*Nv*dGdnp*J;
            end
        end
        
        
        % AUTO INFLUENCE TERM
        % => No Singularity in this case, simply splitting the quadrature
    else
        % For now the number of GP is constant
        % integration over M (not modified)
        for mm=1:N_GAUSS_POINT
            u = GAUSS_POINT(mm);
            Nu = SHAPEFUNC(:,mm)';
            dNu = DSHAPEFUNC(:,mm)';
            Jm = [dNu*x_nodem dNu*y_nodem];
            M = [Nu*x_nodem Nu*y_nodem];
            %
            % The element P is split into two parts with M the common node
            %
            %--------------------------------------------------------------
            % *** First element
            %--------------------------------------------------------------
            N_GAUSS_POINT2 = mm*QuadAutoFactor;
            [GAUSS_POINT2,GAUSS_WEIGHT2] = gauleg(-1,u,N_GAUSS_POINT2);
            
            for pp=1:N_GAUSS_POINT2
                v = GAUSS_POINT2(pp);
                [Nv,dNv] = shape_function(v,Order);
                Jp = [dNv*x_nodep dNv*y_nodep];
                P = [Nv*x_nodep Nv*y_nodep];
                R = norm(P-M);
                %G = -i/4*besselh(0,1,K*R);
                np = [-Jp(2) ; Jp(1)]/norm([-Jp(2) ; Jp(1)]);
                dGdnp = (1i*K/4)*besselh(1,1,K*R)*(P-M)/R*np;
                J = norm(Jm)*norm(Jp);
                % Compute the operator
                Ke = Ke - GAUSS_WEIGHT(mm)*GAUSS_WEIGHT2(pp)*Nu'*Nv*dGdnp*J;
            end
            %--------------------------------------------------------------
            % *** Second element
            %--------------------------------------------------------------
            N_GAUSS_POINT2 = (N_GAUSS_POINT+1 - mm) * QuadAutoFactor;
            [GAUSS_POINT2,GAUSS_WEIGHT2] = gauleg(u,1,N_GAUSS_POINT2);
            
            for pp=1:N_GAUSS_POINT2
                v = GAUSS_POINT2(pp);
                [Nv,dNv] = shape_function(v,Order);
                Jp = [dNv*x_nodep dNv*y_nodep];
                P = [Nv*x_nodep Nv*y_nodep];
                R = norm(P-M);
                %G = -i/4*besselh(0,1,K*R);
                np = [-Jp(2) ; Jp(1)]/norm([-Jp(2) ; Jp(1)]);
                dGdnp = (1i*K/4)*besselh(1,1,K*R)*(P-M)/R*np;
                J = norm(Jm)*norm(Jp);
                % Compute the operator
                Ke = Ke - GAUSS_WEIGHT(mm)*GAUSS_WEIGHT2(pp)*Nu'*Nv*dGdnp*J;
            end
        end
    end
end


%--------------------------------------------------------------------------
%     Single integration (Calling function with a single argument)
%--------------------------------------------------------------------------

if nargin == 1
    for mm=1:N_GAUSS_POINT
        u = GAUSS_POINT(mm);
        Nu = SHAPEFUNC(:,mm)';
        dNu = DSHAPEFUNC(:,mm)';
        J = norm([dNu*x_nodem dNu*y_nodem]);
        M = [Nu*x_nodem Nu*y_nodem];
        
        % Right hand side - incident plane wave
        Phi_inc = exp(1i*K*M(1));
        
        % Compute the operator
        Ke = Ke + 0.5*GAUSS_WEIGHT(mm)*Nu'*Nu*J;
        Fe = Fe + GAUSS_WEIGHT(mm)*Nu'*Phi_inc*J;
    end
end

