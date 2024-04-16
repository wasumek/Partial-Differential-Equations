%% Poisson Problem with 2D FEM (Triangular mesh)
% [Homework 8]   
% Use the Ritz-Galerkin Method with piecewise linear 
% approximation space and nodal basis, to solve 
% the proper variational form of this equation.

%% Report

     clear; read_mesh; bilinear; rhlinear;
     
     % % % % % %       Report     % % % % % % 
    %                                          %
    %            Error is 0.033558             %
    %                                          %
    %   read_mesh: import data file            %
    %   bilinear : assemble stiffness matrix   %
    %   rhlinear : solve r.h.s. equation       %
    %   solver8  : calculate solutions, error, %
    %              and also plot               %
    %                                          %
    %                               by Wasu    %
      % % % % % % % % % % % % % % % % % % % % % 

    
%% calculate the solutions

% an exact solution
exact_sol = zeros(innerV,1);
for i = 1:allV
    exact_sol(i) = f(x(i),y(i));
end

% an appoximated solution
fe_sol = a\rhs';
fe_sol(end+1:numel(x))=0;% 0 at b.c.

%% calculate error(L^2)
sum_err = 0.0;
for alltri = 1:tri_num
    
    v = [ndc(alltri,1) ndc(alltri,2) ndc(alltri,3)];
    
    sum_err = sum_err ... 
              +((exact_sol(v(1))-fe_sol(v(1)))^2 ...
              + (exact_sol(v(2))-fe_sol(v(2)))^2 ...
              + (exact_sol(v(3))-fe_sol(v(3)))^2) ...
              *  vol(alltri)/3;
    
end
err = sqrt(sum_err);
fprintf('\nError is %f\n\n',err);

%% plot with triangle
tri = delaunay(x,y);
trisurf(tri,x,y,fe_sol);
% gtext(0,0,'This is a text you want it to be on the plot');
% waterfall(x,y,fe_sol);
title(sprintf('Poisson Problem (Error=%f)',err))
xlabel('x')
ylabel('y')
zlabel('u_{FE}(x,y)')

