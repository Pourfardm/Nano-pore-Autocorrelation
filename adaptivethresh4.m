function [I2,bw,bw2,label,x,y]=adaptivethresh3(I)
%% What does the function do?
%This funtion do histogram equalization then uses a quadretic exponential...
%smoothness filter to reinforce the circles and then do adaptive thresholding then uses a
%threshold to calculate the centroids
%% The parameters are as below:
%bw is a black and white image containing centroids as white,
%bw2 is the fulfilled pores of the figure (not only the center)
%I2 is an image with red points on the centroids
%label is a labeled picture with bwlabel function so every pore has a distinct label
%x,y are the coordination of the centroid of the figures
%I is the stored original picture
%% 
% clear all;
% load 'H:\MATLAB\nano paper\poori\I.mat' I;
if size(I,3)==3
    I=rgb2gray(I);
end
[n1,n2]=size(I);
% close all;
% imshow(I);       %Original figure

%% 1-histogram equalization
I=histeq(I);
% figure;imshow(I) % Histogram equalized figure

%% 2-Reinforce the circle of the figure
x1 = -2:0.6:2; x2 = x1;
[X1,X2] = meshgrid(x1,x2);
cr=(exp(-(X1.^2+X2.^2)));
k=conv2(single(I),single(cr)/255);%single works with conv2 andhas less storage than double
I3=(uint8(255*k/max(max(k)))); %normalized figure between 0 to 255
I3([1:3,end-2:end],:)=[];
I3(:,[1:3,end-2:end])=[];
im=I3;
[n1,n2]=size(I3);
%% 3-adaptive thresholding: m=mean of 7*7 widows, if I<.7m then black else w
n=7; % n is the size of the smoothness window of adaptive filter
h=im2col(I3,[n,n]);
m=0.8*mean(mean(I3))*ones(size(I3));% 0.8 is a threshod for mean???????
m((n-1)/2+1:n1-(n-1)/2,(n-1)/2+1:n2-(n-1)/2)=reshape(mean(h),n1-n+1,n2-n+1);%m is a smoothed figure
im=(I3<.8*m);%0.8 is empirical because 0.2 of the figure are pores
bw2=im; %bw2 islogical and must be converto to uint8
%figure;imshow(I.*uint8((bw2)))
% figure;imshow(I.*uint8(not(bw2)))


%% 4-finding centroids of remaining areas,elimination of very small
%objects(less than 4 pixels)
label=bwlabel(bw2);
I2(:,:,1)=I;
I2(:,:,2)=I;
I2(:,:,3)=I;
bw=(zeros(size(I))==1);%bw ia logical zero vector
for i=1:max(max(label))
    [x,y]=find(label==i);
    if length(x)>3
    xc=round(mean(x));%centroids of every labeled region( region here are pores)
    yc=round(mean(y));
    I2(xc,yc,1)=255;%we want to read the coordinates of centroids
    I2(xc,yc,2)=0;
    I2(xc,yc,3)=0;
    bw(xc,yc)=1;
    end
end;
% figure;imshow(I2)

%???????????????????????
temp=I;temp(bw2)=255;
%imshow(temp)
I1=uint8(zeros(n1,n2,3));
%I2=I1;
I1(:,:,1)=temp;
I1(:,:,2)=I;
I1(:,:,3)=I;
% figure;imshow(I1); %I1 is a red fulfilled pores
[x,y]=find(bw==1);
end
%% Exececuting the function
% imshow(bw2)
% imshow(bw)
% imshow(I)
% imshow(I2)
% imshow(label2rgb(label))



