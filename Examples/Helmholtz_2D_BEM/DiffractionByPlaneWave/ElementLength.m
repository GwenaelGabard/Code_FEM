function [Le] = ElementLength(Order)

global NODE ELEMENT N_ELEMENT 

N = Order*3;
[point,weight] = gauleg(-1,1,N);
Le = zeros(1,N_ELEMENT);
for k = 1:N_ELEMENT
    x_node = NODE(1,ELEMENT(:,k)).';
    y_node = NODE(2,ELEMENT(:,k)).';
    for m = 1:N
        u = point(m);
        [~,dNu] = shape_function(u,Order);
        Jm = [dNu*x_node dNu*y_node];
        Le(k) = Le(k) + weight(m)*norm(Jm); 
    end
end
