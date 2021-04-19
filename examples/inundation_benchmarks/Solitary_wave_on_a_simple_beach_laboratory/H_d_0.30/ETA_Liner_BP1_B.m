clear all

%--------------------load files--------------

d=load('eta_00165');
n=length(d);
d1=d(1,:);
d=load('eta_00181');
d2=d(1,:);
d=load('eta_00198');
d3=d(1,:);
d=load('eta_00213');
d4=d(1,:);




%--------------------------------------------

dx=0.05;
for i=1:n
    x(i)=60-(i-1)*dx;
end

ii=0;
for i=n:-1:1
    if x(i)>=-10.0
        if x(i)<=20.0
            ii=ii+1;
            kk=i;
            data(ii,2)=d1(i);
            data(ii,3)=d2(i);
            data(ii,4)=d3(i);
            data(ii,5)=d4(i);
        end
    end
end
ll=0;
for i=ii:-1:1
    ll=ll+1;
    data(ll,1)=x(i+kk-1);
end

dlmwrite('NHWAVE_BP4_B.txt', data,'delimiter', '\t','precision', 6);


