
% Open problems
if strcmp(typegeo,'linev')%|strcmp(typegeo,'circlearc') 
    figure; subplot(1,2,1);
    plot(NODE(2,:),abs(U),'r');
    axis square; axis tight;
elseif strcmp(typegeo,'spiral')
    [T,R] = cart2pol(NODE(1,:),NODE(2,:));
    figure; 
    subplot(1,2,1);
    plot(R*2*pi,abs(U),'k');hold on;
    axis tight; axis square;
% Classical      
else
    [T,R] = cart2pol(NODE(1,:),NODE(2,:));      
    T(find(T<0))=2*pi+T(find(T<0)); 
    T=T*180/pi;
    mini = find(T==min(T));
    T = [T( mini:end ) T(1:mini-1) ];
    figure; 
    subplot(1,2,1);
    U = [U( mini:end ); U(1:mini-1)];
    plot(T,real(U),'k');hold on;set(gca,'Xlim',[0 360]);
    Uanal = [Uanal( mini:end ); Uanal(1:mini-1)];
    axis square; axis tight;
    plot(T,real(Uanal),'k.');
    xlabel('Angle: degree');
    title('Analytic numeric comparison');
    axis square; axis tight;
end
      
% Plot the Mesh
subplot(1,2,2);list = [];
for temp=1:N_NODE
    if isempty(find([ELEMENT(1,:) ELEMENT(end,:)]==temp))==0
        list = [list temp];
    end
end
plot(NODE(1,:),NODE(2,:),'k.',NODE(1,list),NODE(2,list),'r*');
title('Geometry');axis equal; axis tight;

