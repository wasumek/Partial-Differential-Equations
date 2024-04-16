
a = sparse(innerV,innerV);
%initialize P
p = [-1 1 0; -1 0 1];

for k = 1:tri_num
    
    v = [ndc(k,1) ndc(k,2) ndc(k,3)];
    
    b11 = y(v(3))-y(v(1));      b12 = -(y(v(2))-y(v(1)));
    b21 = -(x(v(3))-x(v(1)));   b22 = x(v(2))-x(v(1));
    B   = [b11 b12;b21 b22];
    
    for l = 1:3
        i = v(l);
    
        for m = l:3     
            j = v(m);
            
            if (i<=innerV)&&(j<=innerV)
                tmp = ([b11 b12]*p(:,l)*[b11 b12]*p(:,m)...
                      +[b21 b22]*p(:,l)*[b21 b22]*p(:,m))...
                      *(1/(2*det(B)));
               a(i,j) = a(i,j) + tmp;
               a(j,i) = a(i,j);
               
            end
        end
    end
end