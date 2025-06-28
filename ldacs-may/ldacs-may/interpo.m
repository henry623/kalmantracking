function x=interpo(x_in,amp)
x_in=[x_in, x_in(:,1)];
l0=size(x_in,2);
l1=(l0-1)*amp;
x=zeros(size(x_in,1),l1);
for j=1:size(x_in,1)
    x0=x_in(j,:);
    for i=1:l0-1
        x(j,1+(i-1)*amp:i*amp)=(x0(i+1)-x0(i))/amp*[0:amp-1]+x0(i);
    end
end


