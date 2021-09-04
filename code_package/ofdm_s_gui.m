function varargout = ofdm_s_gui(varargin)
% OFDM_S_GUI MATLAB code for ofdm_s_gui.fig
%      OFDM_S_GUI, by itself, creates a new OFDM_S_GUI or raises the existing
%      singleton*.
%
%      H = OFDM_S_GUI returns the handle to a new OFDM_S_GUI or the handle to
%      the existing singleton*.
%
%      OFDM_S_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OFDM_S_GUI.M with the given input arguments.
%
%      OFDM_S_GUI('Property','Value',...) creates a new OFDM_S_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ofdm_s_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ofdm_s_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ofdm_s_gui

% Last Modified by GUIDE v2.5 05-Jan-2020 02:44:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ofdm_s_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @ofdm_s_gui_OutputFcn, ...
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


% --- Executes just before ofdm_s_gui is made visible.
function ofdm_s_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ofdm_s_gui (see VARARGIN)

% Choose default command line output for ofdm_s_gui
handles.output = hObject;
handles.ofdm_real = 0;
handles.inoise = 0;
handles.datanoise_real = 0;
handles.ZZZ_real = 0;
handles.tu2 = 0;
handles.set2 = 0;
handles.set_real_1 = 0;
handles.Z1_1 = 0;
handles.Z1_2 = 0;
handles.Z2_1 = 0;
handles.Z2_2 = 0;
handles.datanoise = 0;
handles.ofdmd = 0;
handles.datanoise3 = 0;
handles.ofdmx = 0;

set(handles.axes2,'visible','off');
axes(handles.axes2);
image_1 = imread('Idea_image.png');
imshow(image_1);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ofdm_s_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ofdm_s_gui_OutputFcn(hObject, eventdata, handles) 
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
%fs=40000;%噪声频率
message_kaishi = 'Simulating... Please Wait！';
icon_2 = 'help';
hhh_1 = msgbox(message_kaishi,'Start the simulation',icon_2);
set(hhh_1,'Position',[700 300 200 80]);% 使用这个语句可以修改msgbox的位置和大小

f1 = str2num(get(handles.edit1,'string'));
fs = str2num(get(handles.edit2,'string'));
c = str2num(get(handles.edit3,'string'));
val1 = get(handles.popupmenu1,'Value');
f_lvboqi = str2num(get(handles.edit5,'string'));

% Baud rate 4000，f1=4e3
f2=64e3; % Standard UART can choose 16 times sampling, can also choose 64 times sampling, I think it should be convenient for frequency division design.
phi=0; jg=f2/f1; jg2=jg/2; a=1; i1=1; i2=1; i3=1; i4=1;
ditong = 2*f_lvboqi/f1; % test best ditong number is 3.05
jieshu = 4; i_mark = 1;% Find an auxiliary pointer to the correct signal
% self-defined head pointer (could use random signal)
baotou = [1 0 0 1 1 1 1 0 0 0 1 1 1 0 1 0 1 0 0 1 1 1 1 0 0 0 1 1 1 0 1 0 1 0 0 1 1 1 1 0 0 0 1 1 1 0 1 0 1 0 0 1 1 1 1 0 0 0 1 1 1 0 1 0];
p_baotou = length(baotou); baotou_cop = (2*baotou-1);

cplex = [1 1; -1 1; -1 -1; 1 -1];% to solve phase error problem
%------------------------------------------------------------------Image processing module
tupian_dizhi = get(handles.edit4,'string');
tu = imread(tupian_dizhi); 
x=rgb2gray(tu);
thresh1 = graythresh(x);
tu2 = im2bw(x,thresh1);
new=(im2double(tu2).*2-1); 
size_1 = size(new);
indata_origin = reshape(new,1,[]);
indata=[baotou indata_origin];% add head pointer
% indata=randsrc(1,2560,[1,-1]);
p=length(indata);
%--------------------------------------------------------------------
%---------------------------------------------------------------------QPSK Modulation Module
% Serial-to-Parallel Conversion
for i=1:1:p
if rem(i,2)==1
    I(i1)=indata(i);
    i1=i1+1;
end
if rem(i,2)==0
    Q(i2)=indata(i);
    i2=i2+1;
end
end

%qpsk modulate
t=1/f2:1/f2:(p*jg2)/f2;
zaibo1=a*cos(2*pi*f1*t-phi);
data1=kron(I,ones(1,jg));
data111=data1.*zaibo1;
zaibo2=a*sin(2*pi*f1*t-phi);
data2=kron(Q,ones(1,jg));
data222=data2.*zaibo2;
s = (data111+data222); 
%-----------------------------------------------------------------------
%------------------------------------------------------------------------OFDM Modulation
%Ifft Prefix
s3 = reshape(s,32,[]);
s4 = ifft(s3,32,1);
s5=[s4(:,end-32+1:end) s4];
s6 = reshape(s5,1,[]);% OFDM Modulated Signal
ofdm_real = s6;

ofdmd=real(s6);
ofdmx=imag(s6);
%--------------------------------------------------------------------------
%---------------------------------------------------------------------------Noise generation and adding module(Directly uncommand one noise to input to the main link)
tt=1/f2:1/f2:((p*jg2)+1024)/f2;
inoise(1,:) = c*sin(2*pi*fs*tt);% Sin Noice
inoise(2,:)=c*square(2*pi*fs*tt,50);% Square Noise (50%)
inoise(3,:)=c*randn(1,length(s6));  % Gaussian Noise
inoise(4,:)=c*sawtooth(2*pi*fs*tt,0.5); % Triangle Noise

inoise(5,:)=0.*ofdmd;% Simplified EFT Noise
for k=1:5:(((p*jg2)+1024))
    inoise(5,k)=c;
end

inoise(6,:)=0.*s6;% Surge Noise
t1=0:1/64e3:158e-6;
for i=1:length(t1)
    if t1(i)<=1.2e-6
        inoise(6,i)=1*exp(0.481*1e12*(t1(i)^2))-1;
         inoise(6,i)=c* inoise(6,i);
    else if t1(i)<130e-6
            inoise(6,i)=1-(t1(i)-1.2e-6)*1e4;
            inoise(6,i)=c*inoise(6,i);
            else if t1(i)>=130e-6
                   inoise(6,i)=-1+(t1(i)-1.2e-6)*1e4-0.570;
                   inoise(6,i)=c*inoise(6,i);
                end
        end
    end
end
for i=11:length(s6)
    inoise(6,i)=inoise(6,i)+inoise(6,(rem(i,11)+1));
end



 datanoise_real = ofdm_real + inoise(val1,:);
 %inoise=inoise+inoise1;
 datanoise=ofdmd+inoise(val1,:);
 datanoise1=0.7*ofdmd+0.3*inoise(val1,:);
 
 datanoise3=ofdmx+inoise(val1,:);
 datanoise4=0.7*ofdmx+0.3*inoise(val1,:);
%--------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------Signal transmission and noise suppression module(ICA）
%Real-Part ICA
s11=datanoise;
s12=datanoise1;

S=[s11;
   s12];
X=S;

Z1=ICA(X);

%datareal=Z(2,:);
% --------------------------------------------------------------------------------
%Imaginary ICA--------------------------------------------------------------------------

s13=datanoise3;
s14=datanoise4;
S=[s13;
   s14];
X=S;
Z2=ICA(X);

handles.Z1_1 = Z1(1,:);
handles.Z1_2 = Z1(2,:);
handles.Z2_1 = Z2(1,:);
handles.Z2_2 = Z2(2,:);
handles.datanoise = datanoise;
handles.ofdmd = ofdmd;
handles.datanoise3 = datanoise3;
handles.ofdmx = ofdmx;

guidata(hObject, handles);

%dataimag=Z(2,:);
% --------------------------------------------------------------------------------Demodulation and sequence correction, phase correction module
for i_zdeofdmd = 1:1:2
    for j_zdeofdmx = 1:1:2
        datareal = Z1(i_zdeofdmd,:);
        dataimag = Z2(j_zdeofdmx,:);
        
        datareturn1=complex(datareal,dataimag);
        datareturn2=complex(-1*datareal,dataimag);
        datareturn3=complex(-1*datareal,-1*dataimag);
        datareturn4=complex(datareal,-1*dataimag);
        dataerror=complex(datanoise,datanoise3);
        
        z_huanyuan_ofdmd(i_zdeofdmd,:) = Z1(i_zdeofdmd,:);
        z_huanyuan_ofdmx(j_zdeofdmx,:) = Z2(j_zdeofdmx,:);
        
        s7 =reshape(datareturn1,32,[]);
        s8=[s7(:,33:end)];
        s9 = fft(s8,32,1);
        s10 = reshape(s9,1,[]);
        
        %Demodulation
        I1=s10.*zaibo1;
        Wc=ditong*f1/f2;                                          %Cut-OFF frequency
        [b,a]=butter(jieshu,Wc,'low');
        I2=filter(b,a,I1);
        
        j=1;
        for iii=jg/2:jg:p*jg/2
            if I2(iii)>=0
                outdataI(j)=1;
            else if I2(iii)<0
                    outdataI(j)=-1;
                end
            end
            j=j+1;
        end
        
        Q1=s10.*zaibo2;    % Q 
        Wc=ditong*f1/f2;                                          %Cut-OFF frequency
        [d,c]=butter(jieshu,Wc,'low');
        Q2=filter(d,c,Q1);
        
        j=1;
        for iii=jg/2:jg:p*jg/2
            if Q2(iii)>=0
                outdataQ(j)=1;
            else if Q2(iii)<0
                    outdataQ(j)=-1;
                end
            end
            j=j+1;
        end
        i3=1;
        i4=1;
        
        for i=1:1:p     
            if rem(i,2)==1
                outdata1(i)=outdataI(i3);
                i3=i3+1;
            end
            if rem(i,2)==0
                outdata1(i)=outdataQ(i4);
                i4=i4+1;
            end
        end
        
   % ef(((fs-fsmin)/step)+1)=error1;
    set1_origin = outdata1;
    baotou_demod(1,:) = set1_origin(1,1:p_baotou);
    
    % ----------------------------------------------------------------------------------
    s7 =reshape(datareturn2,32,[]);
    s8=[s7(:,33:end)];
    s9 = fft(s8,32,1);
    s10 = reshape(s9,1,[]);
    I1=s10.*zaibo1;
    Wc=ditong*f1/f2;                                          
    [b,a]=butter(jieshu,Wc,'low');
    I2=filter(b,a,I1);
    
    j=1;
    for iii=jg/2:jg:p*jg/2
        if I2(iii)>=0
            outdataI(j)=1;
        else if I2(iii)<0
                outdataI(j)=-1;
            end
        end
        j=j+1;
    end
    
    Q1=s10.*zaibo2;    
    Wc=ditong*f1/f2;                                        
    [d,c]=butter(jieshu,Wc,'low');
    Q2=filter(d,c,Q1);
    
    j=1;
    for iii=jg/2:jg:p*jg/2
        if Q2(iii)>=0
            outdataQ(j)=1;
        else if Q2(iii)<0
                outdataQ(j)=-1;
            end
        end
        j=j+1;
    end
    i3=1;
    i4=1;
    
    for i=1:1:p
        if rem(i,2)==1
            outdata3(i)=outdataI(i3);
            i3=i3+1;
        end
        if rem(i,2)==0
            outdata3(i)=outdataQ(i4);
            i4=i4+1;
        end
    end
    
    set3_origin = outdata3;
    baotou_demod(2,:) = set3_origin(1,1:p_baotou);
    
    %-----------------------------------------------------------------------------------
    s7 =reshape(datareturn3,32,[]);
    s8=[s7(:,33:end)];
    s9 = fft(s8,32,1);
    s10 = reshape(s9,1,[]);
    
    I1=s10.*zaibo1;
    Wc=ditong*f1/f2;                                       
    [b,a]=butter(jieshu,Wc,'low');
    I2=filter(b,a,I1);
    
    j=1;
    for iii=jg/2:jg:p*jg/2
        if I2(iii)>=0
            outdataI(j)=1;
        else if I2(iii)<0
                outdataI(j)=-1;
            end
        end
        j=j+1;
    end
    
    Q1=s10.*zaibo2;    
    Wc=ditong*f1/f2;                                         
    [d,c]=butter(jieshu,Wc,'low');
    Q2=filter(d,c,Q1);
   
    j=1;
    for iii=jg/2:jg:p*jg/2
        if Q2(iii)>=0
            outdataQ(j)=1;
        else if Q2(iii)<0
                outdataQ(j)=-1;
            end
        end
        j=j+1;
    end
    i3=1;
    i4=1;
    
    for i=1:1:p
        if rem(i,2)==1
            outdata4(i)=outdataI(i3);
            i3=i3+1;
        end
        if rem(i,2)==0
            outdata4(i)=outdataQ(i4);
            i4=i4+1;
        end
    end
    
    set4_origin = outdata4;
    baotou_demod(3,:) = set4_origin(1,1:p_baotou);
    
    %-----------------------------------------------------------------------------------
    s7 =reshape(datareturn4,32,[]);
    s8=[s7(:,33:end)];
    s9 = fft(s8,32,1);
    s10 = reshape(s9,1,[]);
    
    I1=s10.*zaibo1;
    Wc=ditong*f1/f2;                                         
    [b,a]=butter(jieshu,Wc,'low');
    I2=filter(b,a,I1);
    
    j=1;
    for iii=jg/2:jg:p*jg/2
        if I2(iii)>=0
            outdataI(j)=1;
        else if I2(iii)<0
                outdataI(j)=-1;
            end
        end
        j=j+1;
    end
    
    Q1=s10.*zaibo2;   
    Wc=ditong*f1/f2;                                      
    [d,c]=butter(jieshu,Wc,'low');
    Q2=filter(d,c,Q1);
     
    j=1;
    for iii=jg/2:jg:p*jg/2
        if Q2(iii)>=0
            outdataQ(j)=1;
        else if Q2(iii)<0
                outdataQ(j)=-1;
            end
        end
        j=j+1;
    end
    i3=1;
    i4=1;
    
    for i=1:1:p
        if rem(i,2)==1
            outdata5(i)=outdataI(i3);
            i3=i3+1;
        end
        if rem(i,2)==0
            outdata5(i)=outdataQ(i4);
            i4=i4+1;
        end
    end
    
    set5_origin = outdata5;
    baotou_demod(4,:) = set5_origin(1,1:p_baotou);
    % ----------------------------------------------------------------------------------
    s7 =reshape(dataerror,32,[]);
    s8=[s7(:,33:end)];
    s9 = fft(s8,32,1);
    s10 = reshape(s9,1,[]);

    I1=s10.*zaibo1;
    Wc=ditong*f1/f2;                                          
    [b,a]=butter(jieshu,Wc,'low');
    I2=filter(b,a,I1);
    
    j=1;
    for iii=jg/2:jg:p*jg/2
        if I2(iii)>=0
            outdataI(j)=1;
        else if I2(iii)<0
                outdataI(j)=-1;
            end
        end
        j=j+1;
    end
    
    Q1=s10.*zaibo2;   
    Wc=ditong*f1/f2;                                      
    [d,c]=butter(jieshu,Wc,'low');
    Q2=filter(d,c,Q1);
    
    j=1;
    for iii=jg/2:jg:p*jg/2
        if Q2(iii)>=0
            outdataQ(j)=1;
        else if Q2(iii)<0
                outdataQ(j)=-1;
            end
        end
        j=j+1;
    end
    i3=1;
    i4=1;
    
    for i=1:1:p
        if rem(i,2)==1
            outdata2(i)=outdataI(i3);
            i3=i3+1;
        end
        if rem(i,2)==0
            outdata2(i)=outdataQ(i4);
            i4=i4+1;
        end
    end
    
    errororigin=(sum(abs(indata-outdata2))/length(indata)/2);
    outdata2_real = outdata2(1,p_baotou+1:p);
    set2=reshape(outdata2_real,size_1(1,1),size_1(1,2));
    
    
    %------------------------------------------------------------------------------------Phase Correction Module
    set_zong = [set1_origin;
                set3_origin;
                set4_origin;
                set5_origin;];
    set_zong(1,p+1) = sum(abs(baotou_demod(1,:)-baotou_cop));
    set_zong(2,p+1) = sum(abs(baotou_demod(2,:)-baotou_cop));
    set_zong(3,p+1) = sum(abs(baotou_demod(3,:)-baotou_cop));
    set_zong(4,p+1) = sum(abs(baotou_demod(4,:)-baotou_cop));
    
    set_zong(:,p+2) = i_zdeofdmd;
    set_zong(:,p+3) = j_zdeofdmx;
    
    [xunzhao_v,xunzhao_t] = min(set_zong(:,p+1));
    [xunzhao_k,xunzhao_c] = max(set_zong(:,p+1));
    set_zong(:,p+4) = xunzhao_t;
    set_realorigin_1(i_mark,:) = set_zong(xunzhao_t,:);
        
%     set_realorigin_2 = set_zong(xunzhao_c,:);
%     set_realorigin_2 = -1*set_realorigin_2;
%     error5=(sum(abs(indata-set_realorigin_2))/length(indata)/2);
%     set_realqubaotou_2 = set_realorigin_2(1,p_baotou+1:p);
%     set_real_2 = reshape(set_realqubaotou_2,size_1(1,1),size_1(1,2));
    i_mark = i_mark + 1;
    end
end

 [xunzhao_v,xunzhao_t] = min(set_realorigin_1(:,p+1));
 set_trueorigin = set_realorigin_1(xunzhao_t,:);
 
 i_zdeofdmd_real = set_trueorigin(1,p+2);
 i_zdeofdmx_real = set_trueorigin(1,p+3);
 
 ZZZ_real = complex(cplex(set_trueorigin(1,p+4),1)*z_huanyuan_ofdmd(i_zdeofdmd_real,:), cplex(set_trueorigin(1,p+4),2)*z_huanyuan_ofdmx(i_zdeofdmx_real));
 
 set_realsignal = set_trueorigin(1,1:p);
 errorafter=(sum(abs(indata-set_realsignal))/length(indata)/2);
 set_realqubaotou_1 = set_trueorigin(1,p_baotou+1:p);
 set_real_1 = reshape(set_realqubaotou_1,size_1(1,1),size_1(1,2));
%-------------------------------------------------------------------------------------Plot Utils（Frequency）
set(handles.edit6,'string',100*errororigin);
set(handles.edit7,'string',100*errorafter);

handles.ofdm_real = ofdm_real;
handles.inoise = inoise(val1,:);
handles.datanoise_real = datanoise_real;
handles.ZZZ_real = ZZZ_real;
handles.tu2 = tu2;
handles.set2 = set2;
handles.set_real_1 = set_real_1;
guidata(hObject,handles);

close(hhh_1);
message_jieshu = 'Simulation Over！';
icon_1 = 'help';
hhh = msgbox(message_jieshu,'Simulation Complete',icon_1);
set(hhh,'Position',[700 300 200 80]);% Use this statement to modify the position and size of the MSgbox for the noise spectrum  

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

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure;
subplot(4,1,1);plot(handles.Z1_1);title('Demixing signal 1 (real part)');axis([0 300 -5,5]);
subplot(4,1,2);plot(handles.Z1_2);title('Demixing signal 2 (real part)');axis([0 300 -5,5]);
subplot(4,1,3);plot(handles.datanoise);title('Mixed signal (real part)');axis([0 300 -500,500]);
subplot(4,1,4);plot(handles.ofdmd);title('Original Real Part signal');axis([0 300 -1,1]);
figure;
subplot(4,1,1);plot(handles.Z2_1);title('Demixing signal 1（imaginary part）');axis([0 300 -5,5]);
subplot(4,1,2);plot(handles.Z2_2);title('Demixing signal 2（imaginary part）');axis([0 300 -5,5]);
subplot(4,1,3);plot(handles.datanoise3);title('Mixed signal（imaginary part）');axis([0 300 -500,500]);
subplot(4,1,4);plot(handles.ofdmx);title('Original Imaginary Part（imaginary part）');axis([0 300 -1,1]);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tupian_dizhi = get(handles.edit4,'string');
tu = imread(tupian_dizhi); 
x=rgb2gray(tu);
thresh1 = graythresh(x);
tu2 = im2bw(x,thresh1);

figure;
imshow(tu2);
title('Original Image');

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure;
imshow(handles.set2);
title('Image with noise');

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure;
imshow(handles.set_real_1);
title('De-noise Image');

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 f1 = str2num(get(handles.edit1,'string'));
 fs = str2num(get(handles.edit2,'string'));
 c = str2num(get(handles.edit3,'string'));
 val1 = get(handles.popupmenu1,'Value');
 f2 = 64e3;
 [pinlv,fuzhi]=frequen(handles.ofdm_real,f2); %ofdm Spectrum pic
figure
plot(pinlv,fuzhi);
title('Spectrum diagram of OFDM modulated signal')
xlabel('Frequency/Hz');
ylabel('Amplitude/V/m');

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 f1 = str2num(get(handles.edit1,'string'));
 fs = str2num(get(handles.edit2,'string'));
 c = str2num(get(handles.edit3,'string'));
 val1 = get(handles.popupmenu1,'Value');
 f2 = 64e3;
 [pinlv,fuzhi]=frequen(handles.inoise,f2); 
figure
plot(pinlv,fuzhi);
title('Noise spectrum')
xlabel('Frequency/Hz');
ylabel('Amplitude/V/m');

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 f1 = str2num(get(handles.edit1,'string'));
 fs = str2num(get(handles.edit2,'string'));
 c = str2num(get(handles.edit3,'string'));
 val1 = get(handles.popupmenu1,'Value');
 f2 = 64e3;
 [pinlv,fuzhi]=frequen(handles.datanoise_real,f2); 
figure
plot(pinlv,fuzhi);
title('Signal spectrum diagram after noise superposition')
xlabel('Frequency/Hz');
ylabel('Amplitude/V/m');

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 f1 = str2num(get(handles.edit1,'string'));
 fs = str2num(get(handles.edit2,'string'));
 c = str2num(get(handles.edit3,'string'));
 val1 = get(handles.popupmenu1,'Value');
 f2 = 64e3;
 [pinlv,fuzhi]=frequen(handles.ZZZ_real,f2); 
figure
plot(pinlv,fuzhi);
title('Spectrum diagram of signal after noise suppression')
xlabel('Frequency/Hz');
ylabel('Amplitude/V/m');

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

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 f1 = str2num(get(handles.edit1,'string'));
 fs = str2num(get(handles.edit2,'string'));
 c = str2num(get(handles.edit3,'string'));
 val1 = get(handles.popupmenu1,'Value');

f2=64e3;phi=0;jg=f2/f1;jg2=jg/2;a=1;i1=1;i2=1;i3=1;i4=1;
ditong = 3.05;jieshu = 4;i_mark = 1;
baotou = [1 0 0 1 1 1 1 0 0 0 1 1 1 0 1 0 1 0 0 1 1 1 1 0 0 0 1 1 1 0 1 0 1 0 0 1 1 1 1 0 0 0 1 1 1 0 1 0 1 0 0 1 1 1 1 0 0 0 1 1 1 0 1 0];
p_baotou = length(baotou);
baotou_cop = (2*baotou-1);
cplex = [1 1; -1 1; -1 -1; 1 -1];
tupian_dizhi = get(handles.edit4,'string');
tu = imread(tupian_dizhi); 
x=rgb2gray(tu);
thresh1 = graythresh(x);
tu2 = im2bw(x,thresh1);

new=(im2double(tu2).*2-1); 
size_1 = size(new);

indata_origin = reshape(new,1,[]);
indata=[baotou indata_origin];
p=length(indata);
%--------------------------------------------------------------------
%---------------------------------------------------------------------
for i=1:1:p;
if rem(i,2)==1
    I(i1)=indata(i);
    i1=i1+1;
end
if rem(i,2)==0
    Q(i2)=indata(i);
    i2=i2+1;
end
end
t=1/f2:1/f2:(p*jg2)/f2;
zaibo1=a*cos(2*pi*f1*t-phi);
data1=kron(I,ones(1,jg));
data111=data1.*zaibo1;
zaibo2=a*sin(2*pi*f1*t-phi);
data2=kron(Q,ones(1,jg));
data222=data2.*zaibo2;
s = (data111+data222); 
%-----------------------------------------------------------------------
%------------------------------------------------------------------------
%ifft head pointer
s3 = reshape(s,32,[]);
s4 = ifft(s3,32,1);
s5=[s4(:,end-32+1:end) s4];
s6 = reshape(s5,1,[]);
ofdm_real = s6;
ofdmd=real(s6);
ofdmx=imag(s6);
%--------------------------------------------------------------------------
%---------------------------------------------------------------------------
tt=1/f2:1/f2:((p*jg2)+1024)/f2;
inoise(1,:) = c*sin(2*pi*fs*tt);
inoise(2,:)=c*square(2*pi*fs*tt,50);
inoise(3,:)=c*randn(1,length(s6));  
inoise(4,:)=c*sawtooth(2*pi*fs*tt,0.5);
inoise(5,:)=0.*ofdmd;
for k=1:5:(((p*jg2)+1024))
    inoise(5,k)=c;
end
inoise(6,:)=0.*s6;
t1=0:1/64e3:158e-6;
for i=1:length(t1)
    if t1(i)<=1.2e-6
        inoise(6,i)=1*exp(0.481*1e12*(t1(i)^2))-1;
         inoise(6,i)=c* inoise(6,i);
    else if t1(i)<130e-6
            inoise(6,i)=1-(t1(i)-1.2e-6)*1e4;
            inoise(6,i)=c*inoise(6,i);
            else if t1(i)>=130e-6
                   inoise(6,i)=-1+(t1(i)-1.2e-6)*1e4-0.570;
                   inoise(6,i)=c*inoise(6,i);
                end
        end
    end
end
for i=11:length(s6)
    inoise(6,i)=inoise(6,i)+inoise(6,(rem(i,11)+1));
end

figure
plot(tt,inoise(val1,:));
axis([0 0.001 -c-20 c+20]);
xlabel('Time/s');
ylabel('Field amplitude/V/m');

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



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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
