function varargout = final(varargin)
% FINAL MATLAB code for final.fig
%      FINAL, by itself, creates a new FINAL or raises the existing
%      singleton*.
%
%      H = FINAL returns the handle to a new FINAL or the handle to
%      the existing singleton*.
%
%      FINAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FINAL.M with the given input arguments.
%
%      FINAL('Property','Value',...) creates a new FINAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before final_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to final_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help final

% Last Modified by GUIDE v2.5 08-Jun-2017 11:57:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @final_OpeningFcn, ...
                   'gui_OutputFcn',  @final_OutputFcn, ...
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


% --- Executes just before final is made visible.
function final_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to final (see VARARGIN)

% Choose default command line output for final
handles.output = hObject;
handles.cutconst=0.01;
handles.rL=0.4;
handles.rH=6;%可根据需要效果调整参数
handles.c=4;
handles.d0=8;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes final wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = final_OutputFcn(hObject, eventdata, handles) 
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
%OPENFILE
[filename,pathname]=uigetfile(...
    {'*.bmp;*.jpg;*.png;*.jpeg','image file(*.bmp;*.jpg;*.png;*.jpeg)';...
    '*.*','all file(*.*)'
    }, ...
    'pick an image');
if isequal(filename,0)||isequal(pathname,0)
    return;
end


fpath=[pathname filename];
global S
S=fpath;
imgsrc=imread(fpath);
handles.img=imgsrc;
axes(handles.axes1);
imshow(imgsrc);
guidata(hObject,handles);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes2);
axis off;
cla reset;

axes(handles.axes3);
axis off;
cla reset;

I=handles.img;
I = im2double(I);
fixedSize = 800;
imHeight = size(I,1);
imWidth = size(I,2);
imDiag = sqrt(imWidth^2 + imHeight^2);
downScale = fixedSize / imDiag;
I = imresize(I,downScale,'bicubic');
 I1=I+(handles.cutconst)*I.*(255-I);
J = rgb2gray(I1);

% Detect the 4 dominant lines in the image
BW = edge(J,'canny',[],6);
[H,T,R] = hough(BW);
P = houghpeaks(H,4,'Threshold',0,'NHoodSize',[ceil(fixedSize/4)+1, 41]);
theta = T(P(:,2));
rho = R(P(:,1));
if (length(theta) < 4)
    fprintf('Could not find all lines\n');
    return;
end

% Display detected lines
%figure;
%subplot(1,2,1);
axes(handles.axes2);
axis off;
imagesc(I(:,:,2));
shading flat; hold on; colormap gray; axis equal;
for i = 1:length(theta)
    t = theta(i) / 180 * pi;
    r = rho(i);
    if (abs(t) > pi/4)
        u = 0:size(J,2);
        v = 1 + (r - u*cos(t) )/ sin(t);    
    else
        v = 0:size(J,1);
        u = 1 + (r - sin(t)*v) / cos(t);
    end
    plot(u,v,'g-','LineWidth',2);
end
[~,order] = sort(abs(theta)); 
coefficients = [cos(theta' / 180 * pi), sin(theta' / 180 * pi), -rho'];
L1 = coefficients(order(1),:); % start with smallest angle
L2 = coefficients(order(2),:); % parallel to L1
L3 = coefficients(order(3),:); % perpendicular to L1 and L2
L4 = coefficients(order(4),:); % parallel to L3

% Compute p1, p2, p3, p4 as intersections of L1,...,L4;
% see above figure.
p = zeros(4,3); 
p(1,:) = cross(L1,L3); % p1
p(2,:) = cross(L2,L3); % p2
p(3,:) = cross(L2,L4); % p3
p(4,:) = cross(L1,L4); % p4
p = p(:,1:2) ./ [p(:,3), p(:,3)];

% We now need to reorder the points such that they correspond to 
% the physical document vertices

% Sort points into clockwise order (see above figure)
v = p(2,:) - p(1,:); 
w = p(3,:) - p(1,:);
a = (v(1)*w(2) - w(1)*v(2)) / 2; % signed area of triangle (p1,p2,p3)
if (a < 0)
    p = p(end:-1:1,:); % reverse vertex order
end

% Make sure that first vertex p1 lies either topmost or leftmost,
% and lies at the start of a short side
edgeLen = [
    norm(p(1,:) - p(2,:));
    norm(p(2,:) - p(3,:));
    norm(p(3,:) - p(4,:));
    norm(p(4,:) - p(1,:))
    ];
sortedEdgeLen = sort(edgeLen);
startsShortEdge = (edgeLen' <= sortedEdgeLen(2)); % flag which ones are at a short edge
% find the one that is closest to top left corner of image
idx = find(startsShortEdge,2);
if (norm(p(idx(1),:)) < norm(p(idx(2),:)))
    p1Idx = idx(1);
else
    p1Idx = idx(2);
end
order = mod(p1Idx - 1 : p1Idx + 2, 4) + 1; % reorder to start at right vertex
p = p(order,:);

% Diplay vertex order
for i = 1:4
    text(p(i,1),p(i,2),sprintf('%d',i),'BackgroundColor','white');
end

% Define the document vertices in paper coordinates
% We use the A4 standard: 210x280mm
paperWidth = 210; 
paperHeight = 280;
q = [
    0,          0;               % corresponds to p1
    paperWidth, 0;               % corresponds to p2
    paperWidth, paperHeight;     % corresponds to p3
    0,          paperHeight;     % corresponds to p4
];

% Compute homography p -> q
H = computeHomography(p,q);

% Modify homography to trim the edges by 2% (slight enlargement)
trimPercentage = 2;
eps = 1 / (1 - trimPercentage/100) - 1;
trimEdges = [
    1+eps, 0,    -eps*paperWidth/2;
    0,     1+eps, -eps*paperHeight/2;
    0,     0,     1;
];
H = trimEdges * H;

% Warp original image to paper coordinates
h = size(I,1);
w = size(I,2);
u = repmat(1:w,[h 1]);
v = repmat(1:h,[w 1])';
uv = [u(:),v(:),ones(length(u(:)),1)];
uv = uv*H';
uv = uv(:,1:2) ./ [uv(:,3),uv(:,3)];
u = reshape(uv(:,1),h,w);
v = reshape(uv(:,2),h,w);
%subplot(1,2,2);%output
axes(handles.axes3);
surf(u,v,zeros(h,w),I); 
shading interp; 
colormap gray; 
axis equal;
axis([0 paperWidth 0 paperHeight]);
set(gca,'YDir','reverse');
%
pix=getframe(handles.axes3);
handles.cutimg30=pix.cdata;
guidata(hObject,handles);



function H = computeHomography(p,q)
% Quick implementation of homography computation
% There are better ways to do this, see Hartley & Zisserman's book
A = [p(:,1),p(:,2),ones(4,1),zeros(4,3),-p(:,1).*q(:,1),-p(:,2).*q(:,1);
     zeros(4,3),p(:,1),p(:,2),ones(4,1),-p(:,1).*q(:,2),-p(:,2).*q(:,2)];
b = [q(:,1);q(:,2)];    
Hvec = A\b;
H = reshape([Hvec;1],3,3)';



% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pix=getframe(handles.axes3);

[filename,pathname]=uiputfile({'*.bmp','BMP files';'*.jpg;','JPG files'},...
    'Pick an Image');
if  isequal(filename,0)||isequal(pathname,0)
    return;
else
    fpath=[pathname filename];
    handles.fpath=fpath;
    guidata(hObject,handles);
end

imwrite(pix.cdata,fpath);



% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


guidata(hObject);
pix=handles.cutimg30;
imgsrc1=pix;
imgsrc2=rgb2gray(imgsrc1);
rL=handles.rL;
rH=handles.rH;%可根据需要效果调整参数
c=handles.c;
d0=handles.d0;
imgsrc2=TTLB(imgsrc2,rL,rH,c,d0);
handles.aftttlb=imgsrc2;
guidata(hObject,handles);
%imgsrc2=1.2*imgsrc2;
 %imgsrc3=imgsrc2+0.001*imgsrc2.*(255-imgsrc2);
axes(handles.axes4);
cla reset;
axes(handles.axes4);
imshow(imgsrc2);




% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pix=getframe(handles.axes4);

[filename,pathname]=uiputfile({'*.bmp','BMP files';'*.jpg;','JPG files'},...
    'Pick an Image');
if  isequal(filename,0)||isequal(pathname,0)
    return;
else
    fpath=[pathname filename];
    handles.fpathttlb=fpath;
    guidata(hObject,handles);
end

imwrite(pix.cdata,fpath);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc;
close all;
close(gcf);
clear;

function result=TTLB(I,rL,rH,c,d0)
 
I=double(I);
[M,N]=size(I);

I1=log(I+1);%取对数
FI=fft2(I1);%傅里叶变换
n1=floor(M/2);
n2=floor(N/2);
for i=1:M
    for j=1:N
        D(i,j)=((i-n1).^2+(j-n2).^2);
        H(i,j)=(rH-rL).*(exp(c*(-D(i,j)./(d0^2))))+rL;%高斯同态滤波
    end
end
I2=ifft2(H.*FI);%傅里叶逆变换
I3=real(exp(I2));

I4=I3+0.01*I3.*(255-I3);
bb=imadjust(uint8(I3));

result=bb;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
editstr=get(handles.slider1,'Value');
set(handles.edit1,'String',num2str(editstr));
handles.cutconst=editstr;
guidata(hObject,handles);


% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



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


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
editstr=get(handles.slider2,'Value');
set(handles.edit2,'String',num2str(editstr));
handles.rL=editstr;
guidata(hObject,handles);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
editstr2=get(handles.slider3,'Value');
set(handles.edit3,'String',num2str(editstr2));
handles.rH=editstr2;
guidata(hObject,handles);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
editstr1=get(handles.slider4,'Value');
set(handles.edit4,'String',num2str(editstr1));
handles.c=editstr1;
guidata(hObject,handles);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
editstr=get(handles.slider5,'Value');
set(handles.edit5,'String',num2str(editstr));
handles.d0=editstr;
guidata(hObject,handles);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%imshow(handles.aftttlb);

     businessCard   = handles.aftttlb;
     ocrResults     = ocr(businessCard)
     recognizedText = ocrResults.Text;
        roi = [62 11 900 800];
   axes(handles.axes5);
   cla reset;
    Iocr = insertText(businessCard, roi(1:2), recognizedText, 'AnchorPoint', 'RightTop', 'FontSize',24);
     % Iocr2         = insertObjectAnnotation(Iocr, 'rectangle', ...
           %                ocrResults.WordBoundingBoxes, ...
           %                ocrResults.WordConfidences);
     imshow(Iocr);
     set(handles.edit7,'String',ocrResults.Words);
     %imshow(Iocr2);
     
  
     
    % set(handles.edit6,'String', recognizedText);
     


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s = sprintf('本程序类似Microsoft Office Lens，\n 用于在缺乏专业扫描仪的情况下\n自动分析截取照片并且\n提供图像修正和OCR的功能\n Copyright@ JW. Zhang\n');
msgbox(s,'关于本程序','help')


% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s = sprintf('licensed under the BSD 2-clause "Simplified" License\n Some Codes Refer to Tom Mertens, tom.mertens@gmail.com\n https://github.com/Mericam/document-scanner/ \n Coded by JW Zhang(AmateurZhang)\n https://github.com/AmateurZhang\n');
msgbox(s,'license','help')



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
