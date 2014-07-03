%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEM Method
fprintf('METHOD_BEM : Domain %i (%i elements)\n',ELEMENT_DOMAIN(1,list(1)),length(list));

%--------------------------
% DOUBLE INTEGRATION
%-------------------------
% Loop over the elements
for m=1:length(list)
   % Loop over the elements
   for n=1:length(list)
       % call the elementary routine
       [Ke,Fe,Re,Ve] = feval(char(ELEMENT_NAME(ELEMENT_LIST(list(m)))),list(m),list(n));
       % List of dofs of element m
       dof_listm = DOF_ELEMENT(1:N_DOF_ELEMENT(list(m)),list(m));
       % List of dofs of element n
       dof_listn = DOF_ELEMENT(1:N_DOF_ELEMENT(list(n)),list(n));
       % Assembly
       K(dof_listm,dof_listn) = K(dof_listm,dof_listn) + Ke;
       F(dof_listm) = F(dof_listm)+Fe;
       % Sum up the linear relations
       Le = find(max(abs(Re),[],2));
       R(dof_listm,dof_listn) = R(dof_listm,dof_listn) + Re;
       V(dof_listm) = V(dof_listm) + Ve;
   end
end

%--------------------------
% SIMPLE INTEGRATION
%-------------------------
for m=1:length(list)
   % call the elementary routine
   [Ke,Fe,Re,Ve] = feval(char(ELEMENT_NAME(ELEMENT_LIST(list(m)))),list(m));
   % List of dofs of element m
   dof_list = DOF_ELEMENT(1:N_DOF_ELEMENT(list(m)),list(m));
   % Assembly
   K(dof_list,dof_list) = K(dof_list,dof_list) + Ke;
   F(dof_list) = F(dof_list)+Fe;
   % Sum up the linear relations
   R(dof_list,dof_list) = R(dof_list,dof_list) + Re;
   V(dof_list) = V(dof_list) + Ve;
end

% clear the temp variables
clear Le Fe Ke Re Ve
