clear all ; close all ; clc ;

for k=1:10
    for j=1:15
        for i=1:5
            x(i,j,k)=2*i+10;
            y(i,j,k)=2*j-15;
            z(i,j,k)=2*k-122; 
            v(i,j,k)=i*x(i,j,k)+y(i,j,k)+z(i,j,k); 
        end 
    end 
end
for k=1:12
    for j=1:10
        for i=1:3
            xx(i,j,k)=2*k+10;
            yy(i,j,k)=2*i-15;
            zz(i,j,k)=2*j-11; 
            %v(i,j,k)=i*x(i,j,k)+y(i,j,k)+z(i,j,k); 
        end 
    end 
end
vv=griddata(x,y,z,v,xx,yy,zz);