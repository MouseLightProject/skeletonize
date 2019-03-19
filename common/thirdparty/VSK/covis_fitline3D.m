function [a,b,c,d]=covis_fitline3D(point)
%
% Fit a line in 3D space from a group of skeleton points by using Least Square Method. 
%
% A line in 3D space can be written as following formula:
% X=a*Z+b;      (1)
% Y=c*Z+d;      (2)
%
% Then we can compute the coefficient a,b,c,d by applying LSM to equation(1) 
% and equation(2) respectively.
% 
% input:
% point is a group of skeleton points
%
% output:
% a,b,c,d are coefficient of the line.
%
% ---------------------------
% written by Li Liu in 12/21/2012 
% l.liu6819@gmail.com
%

x=point(:,1);
y=point(:,2);
z=point(:,3);

%transform each line equation into "Aq=b form", and q=inv(A'*A)*A'*b
A=[z,ones(size(z))];
b1=x;
q1=inv(A'*A)*A'*b1;
b2=y;
q2=inv(A'*A)*A'*b2;

a=q1(1);
b=q1(2);
c=q2(1);
d=q2(2);
