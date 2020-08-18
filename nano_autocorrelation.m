% nano autocorrelation code
% 8-2-92
clear;
close all;
clc;
% J = imread('pout.tif');
% J=imread('C:\Users\office\Pictures\Shingubara\Shingu01.tif');
% J=imread('C:\Users\office\Pictures\Shingubara\Shingu02-02-20.bmp');
% J=imread('C:\Users\office\Pictures\SEM-Arab\1-1-2-.tif');imshow(J);
% J=imread('C:\Users\office\Pictures\SEM-Arab\1-1B-2-.tif');imshow(J);
% J=imread('C:\Users\office\Pictures\test test test\3.tif');%imshow(J);1.tif,2.tif
% J=J(:,:,1:3);
% save autocorrelationnano J;
% load autocorrelationnano;
J=imread('special 2.tif');%Perfect
% J=rgb2gray(J);
[I2,bw,bw2,label,x,y]=adaptivethresh4(J);
close all;
I=bw2;
tic
% I=rgb2gray(J);
for i=0:50
    for j=0:50
        A=I;
        B=I;
        A(1:i,:)=[];
        A(:,1:j)=[];
        B(end-i+1:end,:)=[];
        B(:,end-j+1:end)=[];
        II(i+1,j+1)=corr2(A,B);
%         II(i+1,j+1)=mean(mean(A.*B));
    end
    i
end
toc
colormap(jet)
[X,Y] = meshgrid(1:i+1,1:j+1);
JJ=II;
II(1,1)=0;
I_alaki=double(I);
subplot(2,2,1); surf(X,Y,JJ);title('Surf of the autocorrelation matrix');xlabel('x');ylabel('y');zlabel('Autocorrelation');
subplot(2,2,2); imshow(I_alaki);
subplot(2,2,3); surf(X,Y,II);title('Surf of the autocorrelation matrix without (0,0) pixel');xlabel('x');ylabel('y');zlabel('Autocorrelation');
subplot(2,2,4); contour(X,Y,II,100);title('Contour of the autocorrelation matrix without (0,0) pixel');xlabel('x');ylabel('y');
JJ3=JJ; X3=X; Y3=Y; J3=J; I3=I; II3=II;
save autocorrelation3.mat II3 JJ3 X3 Y3 J3 I3;

std2(JJ).^2
%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;
% clc;
close all;
load autocorrelation3;
image=II3;

% image =II;
% image =JJ;
low=0;
high=max(max(image))/2;
a1=circshift(image,[1,0]);
a2=circshift(image,[-1,0]);
a3=circshift(image,[0,1]);
a4=circshift(image,[0,-1]);
a5=circshift(image,[1,1]);
a6=circshift(image,[-1,1]);
a7=circshift(image,[1,-1]);
a8=circshift(image,[-1,-1]);
uuu=max(max(image));
for j=1:10
    low=uuu/10*(j-1);
    high=uuu;
    B=(image>a1)& (image>a2)&(image>a3)&(image>a4)&(image>a5)&(image>a6)&(image>a7)&(image>a8)& (image>=low)& (image<=high);
    imshow(B);
    %     pause;
    number(j)=sum(sum(B==1))
end
% colormap('default')
subplot(2,2,2),imshow(J);
subplot(2,2,1); surf(X,Y,image);
figure, imshow(double(not(I))); 
%%%%%%%%%%%%%%%%
% clear image11;



theta=39;
number11=number;
for i=1:size(image,1)
image11(i)=image(i,round(i*tan(theta*pi/180)));% 60 degree angle
end

numbernew=number;
for i=1:size(image,1)
imagenew(i)=image(i,round(i*tan(theta*pi/180)));% 60 degree angle
end
image11=[image11;imagenew];
number11=[number11;numbernew];

save three.mat image11 number11;
load three;
% number11=[number11;numbernew];
% for g=[1:23]
%     ff(g)=image11(2,g*2)';
% end
ss2=resample(image11(2,:),22,39);
ss3=resample(image11(2,:),22,61);
ss2=[ss2,zeros(1,size(image11,2)-length(ss2))];
ss3=[ss3,zeros(1,size(image11,2)-length(ss3))];
image11(2,:)=ss2;
image11(3,:)=ss3;
% figure,plot(image11(1,:)');
% hold,plot(ff);
% hold,plot(image11(3,:)');


figure,plot(image11');

figure,plot(number11');
save figure11_a image11;
load figure11_a;





% A9	A8	A7	A6	A5	A4	A3	A2	A1	A0	
% 28	13	3	2	1	1	1	1	1	1	a
% 22	15	6	3	3	3	1	1	1	1	b
% 26	24	11	5	3	3	1	1	1	1	c
% 25	18	4	3	3	1	1	1	1	1	d
% 26	15	7	6	2	1	1	1	1	1	e
% 25	12	4	2	1	1	1	1	1	1	f
% 10	9	9	9	7	3	2	1	1	1	1
% 5	4	1	1	1	1	1	1	1	1	2   2
% 10	5	1	1	1	1	1	1	1	1	3
% 20	20	18	12	7	4	1	1	1	1	sem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% fff=[
% 28	13	3	2	1	1	1	1	1	1	
% 22	15	6	3	3	3	1	1	1	1	
% 26	24	11	5	3	3	1	1	1	1	
% 25	18	4	3	3	1	1	1	1	1	
% 26	15	7	6	2	1	1	1	1	1	
% 25	12	4	2	1	1	1	1	1	1	
% 10	9	9	9	7	3	2	1	1	1	
%  5	4	1	1	1	1	1	1	1	1	
% 10	5	1	1	1	1	1	1	1	1	
% 20	20	18	12	7	4	1	1	1	1]
% close all;
% lll=[ fff(3,:); fff(2,:); fff(6,:)]
% plot(lll','--s');
% %%%%%
% close all;
% lll=[fff(7,:); fff(8,:); fff(9,:)]
% plot(lll','--s');