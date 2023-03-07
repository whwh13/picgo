function d=skew_pair(line1,line2)

q1=line1(1,:);
q2=line2(1,:);
t=line1(2,:)-q1;
s=line2(2,:)-q2;
A=[-dot(t,t) dot(s,t);-dot(t,s) dot(s,s)];
B=[dot(q2-q1,t);dot(q2-q1,s)];
%if inverse does exist.
if(abs(det(A))>=1e-10)
    
X=-A\B;
d=norm(q2-q1-t.*X(1)+s.*X(2));
else
    d=0;
end
end