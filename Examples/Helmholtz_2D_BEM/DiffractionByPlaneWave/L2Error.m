function Er = L2Error(U,Ua)

% Computing the L2 error

global ELEMENT NODE N_ELEMENT DOF_ELEMENT
global N_GAUSS_POINT Order

A = 0;
B = 0;

for m = 1:N_ELEMENT
    x_nodem = NODE(1,ELEMENT(:,m))';
    y_nodem = NODE(2,ELEMENT(:,m))';
    [GAUSS_POINTM,GAUSS_WEIGHTM] = gauleg(-1,1,N_GAUSS_POINT);
  
    for mm = 1:N_GAUSS_POINT
        u = GAUSS_POINTM(mm);
        [Nu,dNu] = shape_function(u,Order);
        J = norm([dNu*x_nodem dNu*y_nodem]);
        Phi = Nu*U(DOF_ELEMENT(:,m));
        Phia = Nu*Ua(DOF_ELEMENT(:,m));
        A = A + GAUSS_WEIGHTM(mm)*J*( abs(Phi-Phia)^2  );
        B = B + GAUSS_WEIGHTM(mm)*J*(abs(Phia)^2  );

    end
end

Er = 100*sqrt(A/B);