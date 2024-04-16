rhs = zeros(1,innerV);
vol = zeros(1,innerV);

for i = 1:tri_num
    
    v = [ndc(i,1) ndc(i,2) ndc(i,3)];
    
    b11=y(v(3))-y(v(1));      b12=-(y(v(2))-y(v(1)));
    b21=-(x(v(3))-x(v(1)));   b22=x(v(2))-x(v(1));
    B = [b11 b12;b21 b22];
    
    vol(i) = det(B)/2; % vol. of one triangle
    
    for k = 1:3
        
        if (v(k)<=innerV) 
            rhs(v(k)) = rhs(v(k))...
                        + (vol(i)/3)...
                        * f(x(v(k)),y(v(k)))...
                        * (2.0*pi^2);
        end
            
    end
    
end



