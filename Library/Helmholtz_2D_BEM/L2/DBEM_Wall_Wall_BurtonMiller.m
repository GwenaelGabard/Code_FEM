function [Ke,Fe,Re,Ve] = DBEM_Wall_Wall_BurtonMiller(m,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helmholtz equation BEM formulation (direct BEM)
% Isoparametric high-order element (Lagrange polynomial)
% Interpolation depends on order
% L2->1, L3->2, ...
% Right hand side: Diffraction by a plane wave 
% Formulation:

% Using Burton Miller DBEM formulation with singularity substraction by
% Bonnet
% 1. Burton AJ, Miller GF. The application of integral equation methods to 
% the numerical solution of some exterior boundary-value problems. 1971
% 2. Bonnet M, Guiggiani M. Galerkin BEM with direct evaluation of hyper-
% singular integrals. 2003
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Order

NrOfNodes = (Order+1); 

if nargin==0
   Ke = [1 NrOfNodes ones(1,NrOfNodes) 0];
   return
end

global omega N_GAUSS_POINT alpha
global ELEMENT NODE NODE_DATA
global GAUSS_POINT GAUSS_WEIGHT SHAPEFUNC DSHAPEFUNC 
global N_DOF_ELEMENT

RHO = NODE_DATA(1,ELEMENT(1,m));
C0 = NODE_DATA(2,ELEMENT(1,m));
K = (omega/C0);

Ke = zeros(NrOfNodes,NrOfNodes);
I3S = zeros(NrOfNodes,NrOfNodes);
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
if nargin ==2;
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
            % integrating over P
            M = [Nu*x_nodem Nu*y_nodem]; 
            for pp=1:N_GAUSS_POINT 
                v = GAUSS_POINT(pp);
                Nv = SHAPEFUNC(:,pp)'; 
                dNv = DSHAPEFUNC(:,pp)';
                Jp = [dNv*x_nodep dNv*y_nodep];
                P = [Nv*x_nodep Nv*y_nodep]; 
                R = norm(P-M); 
                nm = [-Jm(2) ; Jm(1)]/norm([-Jm(2) ; Jm(1)]); % Normale sortante du domaine de rayonnement en M
                np = [-Jp(2) ; Jp(1)]/norm([-Jp(2) ; Jp(1)]); 
                G = -1i/4*besselh(0,1,K*R);
                dGdnp = (1i*K/4)*besselh(1,1,K*R)*(P-M)/R*np;              
                J = norm(Jm)*norm(Jp); 
                % Compute the operator
                Ke = Ke + GAUSS_WEIGHT(mm)*GAUSS_WEIGHT(pp)*(- Nu'*Nv*dGdnp*J - alpha*G*( K^2*Nu'*Nv*(nm'*np)*J - dNu'*dNv)); 
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
                % normals
                nm = [-Jm(2) ; Jm(1)]/norm([-Jm(2) ; Jm(1)]); % Normale sortante du domaine de rayonnement en M
                np = [-Jp(2) ; Jp(1)]/norm([-Jp(2) ; Jp(1)]); % Normale sortante du domaine de rayonnement en P
                % kernel
                G = -1i/4*besselh(0,1,K*R);
                dGdnp = (1i*K/4)*besselh(1,1,K*R)*(P-M)/R*np;
                % Jacobian               
                J = norm(Jm)*norm(Jp); 
                % Compute the operator
                f_uv = K^2*Nu'*Nv*(nm'*np)*J - dNu'*dNv;
                % I2
                Ke = Ke + GAUSS_WEIGHT(mm)*GAUSS_WEIGHT2(pp)*( -Nu'*Nv*dGdnp*J);
                % I3
                Ke = Ke - GAUSS_WEIGHT(mm)*GAUSS_WEIGHT2(pp)*alpha*G*f_uv;
                % I3R - regular part -> substracting singularity
				Ke = Ke + GAUSS_WEIGHT(mm)*GAUSS_WEIGHT2(pp)*alpha/(2*pi)*f_uv*log(abs(u-v));
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
                % normals
                nm = [-Jm(2) ; Jm(1)]/norm([-Jm(2) ; Jm(1)]); % Normale sortante du domaine de rayonnement en M
                np = [-Jp(2) ; Jp(1)]/norm([-Jp(2) ; Jp(1)]); % Normale sortante du domaine de rayonnement en P
                % Kernel
                G = -1i/4*besselh(0,1,K*R);
                dGdnp = (1i*K/4)*besselh(1,1,K*R)*(P-M)/R*np;
                % Jacobian               
                J = norm(Jm)*norm(Jp); 
                % Compute the operator
                f_uv = K^2*Nu'*Nv*(nm'*np)*J - dNu'*dNv;
                % I2
                Ke = Ke + GAUSS_WEIGHT(mm)*GAUSS_WEIGHT2(pp)*( -Nu'*Nv*dGdnp*J);
                % I3
                Ke = Ke - GAUSS_WEIGHT(mm)*GAUSS_WEIGHT2(pp)*alpha*G*f_uv;
                % I3R - regular part -> substracting singularity
				Ke = Ke + GAUSS_WEIGHT(mm)*GAUSS_WEIGHT2(pp)*alpha/(2*pi)*f_uv*log(abs(u-v));
            end
        end
        
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Computing the singular part
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%-alpha/(2*pi)* log(u-v)*f(u,v)
		n_gauss_point_KSI = N_GAUSS_POINT*QuadAutoFactor;
		n_gauss_point_ZETA = n_gauss_point_KSI;
		
		% Classical GP and Specific Log GP from 0 to 1
		[nodesNor,weightsNor] = gauleg(0,1,12);
		nodesLog  = [ 0.006548722279080035, 0.03894680956045022, 0.0981502631060046, 0.1811385815906331, 0.2832200676673157, 0.398434435164983, 0.5199526267791299, 0.6405109167754819, 0.7528650118926111, 0.850240024421055, 0.926749682988251, 0.977756129778486];
		weightsLog = [-0.09319269144393, -0.1497518275763289, -0.166557454364573, -0.1596335594369941, -0.1384248318647479, -0.1100165706360573, -0.07996182177673273, -0.0524069547809709, -0.03007108900074863, -0.01424924540252916, -0.004899924710875609, -0.000834029009809656];
		
		epsi = (24/N_GAUSS_POINT)/2;
		        
		% ksi : -1 -> 1
		[gauss_point_KSI,gauss_weight_KSI] = gauleg(-1,1,n_gauss_point_KSI);
		% zetaa goes from 0 to 2 (in 2 steps)
        [gauss_point_ZETA,gauss_weight_ZETA] = gauleg(epsi,2,n_gauss_point_ZETA);
        
        % integration over ksi from 0 to 2
		for mm=1:n_gauss_point_KSI		
            ksi = gauss_point_KSI(mm) ;
            %-----------------------------------------------------
			% LOG GAUSS POINTS (change of variable, we integrate from 0 to 1)
            % zeta: 0 -> epsi
            for pp=1:12
                I3S = I3S + gauss_weight_KSI(mm)*epsi*( fobonnet(epsi*nodesLog(pp),ksi,m)*weightsLog(pp) + fobonnet(epsi*nodesNor(pp),ksi,m)*log(epsi)*weightsNor(pp)  );
            end
            %-----------------------------------------------------
			% SECOND PART, CLASSICAL GAUSS QUADRATURE
			% zeta: epsi -> 2-ksi
            for pp=1:n_gauss_point_ZETA
                zetaa = gauss_point_ZETA(pp);
                I3S = I3S + gauss_weight_KSI(mm)*gauss_weight_ZETA(pp)*log(zetaa)*fobonnet(zetaa,ksi,m);
            end
		end
		
		% adding the term -int(int(alpha/(2.d0*pi)*f_uv*dlog(dabs(u-v)),u,-1,1),v,-1,1)
		Ke = Ke - alpha/(2*pi)*I3S;
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
        J = ([dNu*x_nodem dNu*y_nodem]);
        nm = [-J(2) ; J(1)]/norm([-J(2) ; J(1)]); % Normale sortante du domaine de rayonnement en M
        J = norm(J);        
        M = [Nu*x_nodem Nu*y_nodem];
        
        % Right hand side - incident plane wave
        Phi_inc = exp(1i*K*M(1));
        dPhi_incdnm = nm.'*[1i*K*Phi_inc ; 0];
                
        % Compute the operator
        Ke = Ke + 0.5*GAUSS_WEIGHT(mm)*(Nu'*Nu)*J;
        Fe = Fe + GAUSS_WEIGHT(mm)*Nu'*(Phi_inc + alpha*dPhi_incdnm)*J;
    end
end 
    
end
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
function [out] = fobonnet(zetaa,ksi,m)

global omega Order NODE ELEMENT 

x_nodem = NODE(1,ELEMENT(:,m))';
y_nodem = NODE(2,ELEMENT(:,m))';

K = (omega/340);

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%! Computing f(u(ksi,zetaa),v(ksi,zetaa))+ f(v(ksi,zetaa),u(ksi,zetaa))
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
u = ksi - zetaa/2 - ksi*zetaa/2;
v = ksi + zetaa/2 - ksi*zetaa/2;

% Geometry
[Nu,dNu] = shape_function(u,Order);
[Nv,dNv] = shape_function(v,Order);
Jm = [dNu*x_nodem dNu*y_nodem];  
Jp = [dNv*x_nodem dNv*y_nodem]; 

% Normals
nm = [-Jm(2) ; Jm(1)]/norm([-Jm(2) ; Jm(1)]); 
np = [-Jp(2) ; Jp(1)]/norm([-Jp(2) ; Jp(1)]); 

% Points
M = [Nu*x_nodem Nu*y_nodem];
P = [Nv*x_nodem Nv*y_nodem];

% Jacobien
J = norm(Jm)*norm(Jp); 
    
fuv = K^2*Nu'*Nv*(nm'*np)*J - dNu'*dNv;
fvu = K^2*Nv'*Nu*(nm'*np)*J - dNv'*dNu;

out = (fuv + fvu)*(1-zetaa/2); 

end  











