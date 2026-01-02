system(date)
add(R,S,A)={local(id,x1,y1,x2,y2,x3,y3,lam);id=[0];if(R==id,return(S));if(S==id,return(R));x1=R[1];y1=R[2];x2=S[1];y2=S[2];if(x1==x2,if(y1==-y2,return(id));if(y1==y2,if(y1==0,return(id));lam=(3*x1^2+A)/(2*y1););,lam=(y2-y1)/(x2-x1););x3=lam^2-x1-x2;y3=lam*(x1-x3)-y1;return([x3,y3])}
neg(R)={local(id);id=[0];if(R==id,id,[R[1],-R[2]])}
mP(A,B,m,P)={local(id,res,Q,k);id=[0];res=id;Q=P;k=m;if(k<0,k=-k;Q=neg(Q));while(k>0,if(k%2==1,res=add(res,Q,A));Q=add(Q,Q,A);k=k\2;);res}
q2mat(q,z)={local(i,j,Q);Q=matdiagonal(vector(#z,i,polcoeff(q,2,z[i])));for(i=1,#z-1,for(j=i+1,#z,Q[i,j]=polcoeff(polcoeff(q,1,z[i]),1,z[j])/2;Q[j,i]=Q[i,j];););Q}
glnz(n,B)={local(g,i,r1,r2,c,t);g=matid(n);for(i=1,B,r1=random(n)+1;r2=random(n)+1;c=random(2)*2-1;for(i=1,n,t=g[r1,i];g[r1,i]+=c*g[r2,i];g[r2,i]=t;););g}


#
default(parisizemax,16*1024^3);
default(threadsizemax,16*1024^3);
n=1176817729
m=41

m=11
n=11171^3
e1=ellinit([0,0,0,-m*n,0])
p=5
f=elldivpol(e1,p);
ff=factor(f);
f0=polcyclo(p);
f1=ff[1,1];
f2=ff[#ff[,1],1];
f3=polresultant(f2,y^2-x,x);
f4=polcompositum(f3,y^2+1)[1];
yi=nfgaloisconj(f4);

x1=nfisincl(f2,f4)[1];
y1=lift(polcoeff(factor(Mod(x^2-(x1^3-m*n*x1),f4))[1,1],0,x));
lift(Mod(y1^2-(x1^3-m*n*x1),f4))
P1=[x1,y1];
x2=-x1;
y2=lift(polcoeff(factor(Mod(x^2-(x2^3-m*n*x2),f4))[1,1],0,x));
P2=[x2,y2];
lift(Mod(y2^2-(x2^3-m*n*x2),f4))

quad_pairs=[];
for(i=1,p,for(j=i,p,quad_pairs=concat([quad_pairs,[[i,j]]]);););
quad_pairs


points = [];
system(date)
m1P = vector(p,i,mP(-m*n,0,i-1,Mod([x1,y1],f4)));
n1P = vector(p,i,mP(-m*n,0,i-1,Mod([x2,y2],f4)));
for(m1=1,p,for(n1=1,p,P3=add(m1P[m1],n1P[n1],A);if(P3!=0,P3=lift(P3));points=concat(points,[P3]);print1(".");););
print();
points_affine = points[2..#points];

bb = [1,x,y,x^2,x*y];
[r,c]=[#points_affine, #quad_pairs];
M = matrix(r,c);
for(r1=1,r,Q=points_affine[r1];[x1,y1]=Q;for(c1=1,c,[i,j]=quad_pairs[c1];v=substvec(bb[i],['x,'y],[Mod(x1,f4),Mod(y1,f4)]) * substvec(bb[j],['x,'y],[Mod(x1,f4),Mod(y1,f4)]);M[r1,c1]=lift(v);print1(".");););
print();
system(date)
k1=lift(matker(Mod(M,f4)))
d=matsize(k1);d

z=vector(p,i,eval(Str("z"i)));

qq = [];
for(k=1,d[2],q1=0;idx=1;v=k1[,k];for(i=1,p,for(j=i,p,if(v[idx]!=0,q1+=v[idx]*z[i]*z[j];);idx+=1;););qq=concat(qq,q1););

q2mat(2*qq[1],z)
q2mat(2*qq[2],z)
q2mat(2*qq[3],z)
q2mat(2*qq[4],z)
q2mat(2*qq[5],z)

system(date);


