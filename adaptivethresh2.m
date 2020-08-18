function [I2,bw,bw2,xnbr,ynbr,dist,label,x,y]=adaptivethresh2(I,nn,mm)
%bw is a black and white image containing centroids as white,
%I2 is an image with red points on the centroids,
%nn,mm= thresholds of dist, 7 ,13 for initial images
if size(I,3)==3
    I=rgb2gray(I);
end
[n1,n2]=size(I);
% close all

%1-histogram equalization
I=histeq(I);
figure;imshow(I)

%2-adaptive thresholding: m=mean of 7*7 widows, if I<.7m then black else w
x1 = -2:0.6:2; x2 = x1;
[X1,X2] = meshgrid(x1,x2);
cr=(exp(-(X1.^2+X2.^2)))
k=conv2(single(I),single(cr)/255);
I3=(uint8(255*k/max(max(k))));
I3([1:3,end-2:end],:)=[];
I3(:,[1:3,end-2:end])=[];
im=I3;
[n1,n2]=size(I3);
n=7;
h=im2col(I3,[n,n]);
m=0.8*mean(mean(I3))*ones(size(I3));
m((n-1)/2+1:n1-n+1+(n-1)/2,(n-1)/2+1:n2-n+1+(n-1)/2)=reshape(mean(h),n1-n+1,n2-n+1);
im=(I3<.8*m);
bw2=im;
figure;imshow(I.*uint8(not(bw2)))


%4-finding centroids of remaining areas,elimination of very small
%objects(less than 4 pixels)
label=bwlabel(bw2);
I2(:,:,1)=I;
I2(:,:,2)=I;
I2(:,:,3)=I;
bw=(zeros(size(I))==1);
for i=1:max(max(label))
    [x,y]=find(label==i);
    if length(x)>3
    xc=round(mean(x));
    yc=round(mean(y));
    I2(xc,yc,1)=255;
    I2(xc,yc,2)=0;
    I2(xc,yc,3)=0;
    bw(xc,yc)=1;
    end
end;
figure;imshow(I2)

%
temp=I;temp(bw2)=255;
I1=uint8(zeros(n1,n2,3));
I2=I1;
I1(:,:,1)=temp;
I1(:,:,2)=I;
I1(:,:,3)=I;
figure;imshow(I1);

%5- finding neighbors in a k*k window
k=13;
[x,y]=find(bw==1);
xnbr=zeros(20,length(x));
ynbr=zeros(20,length(x));
dist=zeros(20,length(x));
for i=1:length(x)
    mask=zeros(n1,n2);
    mask(max(x(i)-k,1):min(x(i)+k,n1),max(y(i)-k,1):min(y(i)+k,n2))=1;
    mask(x(i),y(i))=0;
    mask=(mask==1);
    [xn,yn]=find(bw.*mask==1);
    xnbr(1:length(xn),i)=xn;
    ynbr(1:length(xn),i)=yn;
    dist(1:length(xn),i)=sqrt((x(i)-xn).^2+(y(i)-yn).^2);
end
dd=dist(:);dd(dd==0)=[];
figure;hist(dd,30); 
title('deside about thresholds on d')
%hist(sum(not(dist==0)),13);
[x,y]=find(bw==1);
xnbr((dist<nn)|(dist>mm))=0;
ynbr((dist<nn)|(dist>mm))=0;
dist((dist<nn)|(dist>mm))=0;
for i=1:length(x)
    xn=xnbr(:,i);
    yn=ynbr(:,i);
    dn=dist(:,i);
    xn(xn<=0.1)=[];
    yn(yn<=0.1)=[];
    dn(dn<=0.1)=[];
    n=length(xn);
    xnbr(1:n,i)=xn;
    ynbr(1:n,i)=yn;
    dist(1:n,i)=dn;
    xnbr(n+1:end,i)=0;
    ynbr(n+1:end,i)=0;
    dist(n+1:end,i)=0;
end


