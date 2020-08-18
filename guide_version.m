function varargout = guide_version(varargin)
% GUIDE_VERSION MATLAB code for guide_version.fig
%      GUIDE_VERSION, by itself, creates a new GUIDE_VERSION or raises the existing
%      singleton*.
%
%      H = GUIDE_VERSION returns the handle to a new GUIDE_VERSION or the handle to
%      the existing singleton*.
%
%      GUIDE_VERSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIDE_VERSION.M with the given input arguments.
%
%      GUIDE_VERSION('Property','Value',...) creates a new GUIDE_VERSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guide_version_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guide_version_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guide_version

% Last Modified by GUIDE v2.5 27-Jun-2016 15:31:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guide_version_OpeningFcn, ...
                   'gui_OutputFcn',  @guide_version_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before guide_version is made visible.
function guide_version_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guide_version (see VARARGIN)

% Choose default command line output for guide_version
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guide_version wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = guide_version_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.*';},'mytitle','test\');
handles.address=strcat(pathname,filename);
J=imread(handles.address);
axes(handles.axes1);
imshow(J);
save nano.mat J;
guidata(hObject, handles);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load nano.mat;
[I2,bw,bw2,label,x,y]=adaptivethresh4(J);
I=bw2;
axes(handles.axes2);
imshow(I);title('Thresholded version');
drawnow;
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
colormap(jet);
[X,Y] = meshgrid(1:i+1,1:j+1);
JJ=II;
II(1,1)=0;
I_alaki=double(I);
axes(handles.axes3);
surf(X,Y,II);title('Surf of the autocorrelation matrix');xlabel('x');ylabel('y');zlabel('Autocorrelation');
axes(handles.axes4);
contour(X,Y,II,100);title('Contour of the autocorrelation');xlabel('x');ylabel('y');drawnow;
% JJ3=JJ; X3=X; Y3=Y; J3=J; I3=I; II3=II;
set(handles.edit1,'string',num2str(std2(II)));
%%%%%%%%%%%%%%%%%

JJ3=JJ; X3=X; Y3=Y; J3=J; I3=I; II3=II;
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
%     imshow(B);
    %     pause;
    number(j)=sum(sum(B==1))
end
theta=39;
number11=number;
for i=1:size(image,1)
image11(i)=image(i,round(i*tan(theta*pi/180)));% 60 degree angle
end
% pause;
figure,plot(image11');title('A sample image from theta=39 degree');

figure,plot(number11');title('Number of peaks above ten slots of ACF');


guidata(hObject, handles);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load nano.mat;
I=imresize(J,0.5);
clear J;
% load I
template_size=17;%17 is the window size of the template
f=@(x,y,xc,yc,sig)(exp(-((x-xc).^2+(y-yc).^2)/sig)); %function handle
x=-template_size:template_size;y=x;
[x,y]=meshgrid(x,y);
threshold=32;%Why 32 ??????????lower bound of correlation
r=5; %by changing r we have different domain size
for theta0=0:59
    s=f(x,y,0,0,5); % Sigma
    for i=1:6
        theta=i*pi/3;
        xc=r*cos(theta+theta0*pi/180);yc=r*sin(theta+theta0*pi/180);
        s=s+f(x,y,xc,yc,5);
    end
%         imshow(uint8(255*s));xlabel('6 lobes Gaussian');  %For rotating the template
%         [XXX, YYY]=meshgrid(1:size(s,1),1:size(s,2));
%         mesh(XXX,YYY,s); colormap('jet'); pause
        
    J(:,:,theta0+1)=xcorr2(I(:,:,1),s)/(template_size^2);
end
% figure;
[m_val,k]=max(J,[],3);%maximum on the 3rd dimension
% imshow(uint8(8*min(J,[],3)))
% figure;
% imshow(m>32)
mm=200;
% I=rgb2gray(I);
Kr=mm*trimf(k/60,[0,.25,1])+256-mm;%the same size (x,y axis) of J
Kg=mm*trimf(k/60,[0,.75,1])+256-mm;
Kb=mm*max(trimf(k/60,[-1,0,.5]),...
    trimf(k/60,[.5,1,2]))+256-mm;
I=padarray(I,[(template_size),(template_size)]);%to adapt with Kr
I5=I;
II1=uint8(double(I5(:,:,1)).*Kr/255);
II2=uint8(double(I5(:,:,2)).*Kg/255);
II3=uint8(double(I5(:,:,3)).*Kb/255);
I5(:,:,1)=II1;I5(:,:,2)=II2;I5(:,:,3)=II3;
I5(1:template_size,:,:)=[];I5(end-template_size+1:end,:,:)=[];
I5(:,1:template_size,:)=[];I5(:,end-template_size+1:end,:)=[];
axes(handles.axes2);
imshow(I5);title(['Without threshold and r=',num2str(r)]);
I4=I;
II1=I4(:,:,1);
II1(m_val<threshold)=Kr(m_val<threshold);%threshold is lower bound of the corelation
II2=I4(:,:,2);
II2(m_val<threshold)=Kg(m_val<threshold);
II3=I4(:,:,3);
II3(m_val<threshold)=Kb(m_val<threshold);
true_index=m_val<threshold;%index of non-ordered areas
k(true_index)=61;%Assign 61 to non-ordered area
I4(:,:,1)=II1;I4(:,:,2)=II2;I4(:,:,3)=II3;
I4(1:template_size,:,:)=[];I4(end-template_size+1:end,:,:)=[];
I4(:,1:template_size,:)=[];I4(:,end-template_size+1:end,:)=[];
% axes(handles.axes3);
% imshow(I4);title('different domains orientation with different colors with threshold');%title('r=9');
% pause(1);
%% unsupervised segmentqtion
k(1:template_size,:,:)=[];k(end-template_size+1:end,:,:)=[];%k is index of maximum correlation
k(:,1:template_size,:)=[];k(:,end-template_size+1:end,:)=[];
I(1:template_size,:,:)=[];I(end-template_size+1:end,:,:)=[];
I(:,1:template_size,:)=[];I(:,end-template_size+1:end,:)=[];
for i=0:59
    KK(i+1)=numel(find(k==i));
end
figure,plot(KK);title('Original bin');

bin=KK;
n1=5;
shift1=3;
n2=5;
shift2=3;

filt1=hamming(n1);filt1=filt1/sum(filt1);
bin1=conv(bin,filt1);
bin1=[bin1(shift1:end),bin1(1:shift1-1)];
%     plot(bin1,'r');
% n2=20;


filt2=hamming(n2);filt2=filt2/sum(filt2);
bin2=conv(bin1,filt2);
bin2=[bin2(shift2:end),bin2(1:shift2-1)];
hold on;
figure,
plot(bin2,'r');
bias=1;
bin2min=[bin2(end),bin2(1:end-bias)];% one pixel to the right
bin2plus=[bin2(bias+1:end),bin2(bias)];% one pixel to the left
bias=2;
bin2min2=[bin2(end-bias+1:end),bin2(1:end-bias)];% two pixel to the right
bin2plus2=[bin2(bias+1:end),bin2(1:bias)];% two pixel to the left
bias=3;
bin2min3=[bin2(end-bias+1:end),bin2(1:end-bias)];
bin2plus3=[bin2(bias+1:end),bin2(1:bias)];
bias=4;
bin2min4=[bin2(end-bias+1:end),bin2(1:end-bias)];
bin2plus4=[bin2(bias+1:end),bin2(1:bias)];
binfinal=( (bin2>0)...
    &(bin2>=bin2min)   & (bin2>=bin2plus)  & (bin2>=bin2min2) & (bin2>=bin2plus2)...
    & (bin2>=bin2min3) & (bin2>=bin2plus3) & (bin2>=bin2min4) & (bin2>=bin2plus4));
%         & (bin2>=bin2min5) & (bin2>=bin2plus5) & (bin2>=bin2min6) & (bin2>=bin2plus6)...
%         & (bin2>=bin2min7) & (bin2>=bin2plus7) & (bin2>=bin2min8) & (bin2>=bin2plus8));

figure,plot(find(binfinal==1),bin2(binfinal),'mo');
candidate=find(binfinal==1);

%% new label
for iter=1:60
    %     [v INDEX]=min((iter-ind).^2);
    [v INDEX]=min((iter-candidate).^2);
    new_label(iter)=INDEX; %assign each label to its nearest neighbor
end


m=1;
zaviyeh=numel(candidate);
new_label(61)=zaviyeh+1;%non-ordered areas
zaviyeh=zaviyeh+1;

Krr=m*trimf(0:zaviyeh-1,[0,.25,1]*zaviyeh);
Kgg=m*trimf(0:zaviyeh-1,[0,.75,1]*zaviyeh);
Kbb=m*max(trimf(0:zaviyeh-1,[-1,0,.5]*zaviyeh),...
    trimf(0:zaviyeh-1,[.5,1,2]*zaviyeh));% for two part trimf


I6=I;
II1=uint8(double(I6(:,:,1)).*Krr(new_label(k)));
II2=uint8(double(I6(:,:,2)).*Kgg(new_label(k)));
II3=uint8(double(I6(:,:,3)).*Kbb(new_label(k)));
I6(:,:,1)=II1;I6(:,:,2)=II2;I6(:,:,3)=II3;
I6(1:template_size,:,:)=[];I6(end-template_size+1:end,:,:)=[];
I6(:,1:template_size,:)=[];I6(:,end-template_size+1:end,:)=[];
num_domain_reduced=numel(candidate)+1;
axes(handles.axes3);
imshow(I6);title(['unsupervised segmentation with ',num2str(num_domain_reduced),' domains']);

%% Merge regions of image
k_new=new_label(k);
xt1=num_domain_reduced+1;%number of domains 
xt2=0;%number of domains after merging

domain_num=0;
final_domain=zeros(size(I,1),size(I,2));

for domain=1:xt1
    I_color_test=(k_new==domain);
    % labeling
    [L, num] = bwlabel(I_color_test,4);
    final_domain(L~=0)=final_domain(L~=0)+L(L~=0)+domain_num;
    domain_num=domain_num+num;
%      
end
xt2=domain_num;

m=1;
zaviyeh=domain_num;

Krrr=m*trimf(0:zaviyeh-1,[0,.25,1]*zaviyeh);
Kggg=m*trimf(0:zaviyeh-1,[0,.75,1]*zaviyeh);
Kbbb=m*max(trimf(0:zaviyeh-1,[-1,0,.5]*zaviyeh),...
    trimf(0:zaviyeh-1,[.5,1,2]*zaviyeh));% for two part trimf

I_final(:,:,1)=double(I(:,:,1)).*Krrr(final_domain);
I_final(:,:,2)=double(I(:,:,2)).*Kggg(final_domain);
I_final(:,:,3)=double(I(:,:,3)).*Kbbb(final_domain);
I_final=uint8(I_final);

axes(handles.axes4);
imshow(I_final),title(['Final domain with ',num2str(domain_num),' domains']);

%% largest domain size
largest_domain_size=0;
for i=1:domain_num
    temp=numel(find(final_domain==i));
if(temp>largest_domain_size)
    largest_domain_size=temp;
    best_i=i;
end
end
largest_domain=(final_domain==best_i);
figure,
imshow(largest_domain*255);title(['largest domain size with area of ',num2str(largest_domain_size/numel(final_domain))]);




guidata(hObject, handles);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load nano.mat;
I=J;
clear J;
% pause;
[I2,bw,bw2,xnbr,ynbr,dist,label,x,y]=adaptivethresh2(I,7,11);
% close all;
axes(handles.axes2);
imshow(label);drawnow;
[xtri,ytri,xc,yc,trireg,tridevind]=triangulation_ab(I,x,y,xnbr,ynbr,0.2,0.2);
% % previous approach:
% [ntrilabel,mrgim,labelcolim,labelim]=spreading(I,xc,yc,xtri,ytri,trireg);

%my algorithm:
[theta,dtheta,mrgim,label,indv,devind,labelcolim,labelim,corestab]=...
    myangle4(I,bw,x,y,xnbr,ynbr,dist,trireg,tridevind,xtri,ytri,0.001);

[nlabel,mrgim,labelcolim,labelim]=myspreading(I,label,x,y,xnbr,ynbr);
axes(handles.axes3);
imshow(mrgim);

%color coding:
tehta(theta==0)=1;
theta(theta==100)=0;
imerg(I,theta,x,y,xnbr,ynbr);
%bar:
a=[1:59];a=repmat(a,20,1);
Kr=200*trimf(a/60,[0,.25,1])+256-200;
Kg=200*trimf(a/60,[0,.75,1])+256-200;
Kb=200*max(trimf(a/60,[-1,0,.5]),trimf(a/60,[.5,1,2]))+256-200;
J(:,:,1)=Kr;J(:,:,2)=Kg;J(:,:,3)=Kb;
axes(handles.axes4);
imshow(uint8(J))

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
