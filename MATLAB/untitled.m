function varargout = untitled(varargin)
% UNTITLED M-file for untitled.fig
%      UNTITLED, by itself, creates a new UNTITLED or raises the existing
%      singleton*.
%
%      H = UNTITLED returns the handle to a new UNTITLED or the handle to
%      the existing singleton*.
%
%      UNTITLED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNTITLED.M with the given input arguments.
%
%      UNTITLED('Property','Value',...) creates a new UNTITLED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help untitled

% Last Modified by GUIDE v2.5 06-Sep-2012 18:59:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @untitled_OpeningFcn, ...
    'gui_OutputFcn',  @untitled_OutputFcn, ...
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

% --- Executes just before untitled is made visible.
function untitled_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled (see VARARGIN)

% Choose default command line output for untitled
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using untitled.
if strcmp(get(hObject,'Visible'),'off')
    
end

% UIWAIT makes untitled wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = untitled_OutputFcn(hObject, eventdata, handles)
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

cla;
data_set_1x = getappdata(0,'dataset_1x');
data_set_1y = getappdata(0,'dataset_1y');
data_set_1d = getappdata(0,'dataset_1d');
data_set_2x = getappdata(0,'dataset_2x');
data_set_2y = getappdata(0,'dataset_2y');
data_set_2d = getappdata(0,'dataset_2d');
data_set_3x = getappdata(0,'dataset_3x');
data_set_3y = getappdata(0,'dataset_3y');
data_set_3d = getappdata(0,'dataset_3d');
data_set_4x = getappdata(0,'dataset_4x');
data_set_4y = getappdata(0,'dataset_4y');
data_set_4d = getappdata(0,'dataset_4d');

c_data_set_1x = getappdata(0,'c_dataset_1x');
c_data_set_1y = getappdata(0,'c_dataset_1y');
c_data_set_1d = getappdata(0,'c_dataset_1d');
c_data_set_2x = getappdata(0,'c_dataset_2x');
c_data_set_2y = getappdata(0,'c_dataset_2y');
c_data_set_2d = getappdata(0,'c_dataset_2d');
c_data_set_3x = getappdata(0,'c_dataset_3x');
c_data_set_3y = getappdata(0,'c_dataset_3y');
c_data_set_3d = getappdata(0,'c_dataset_3d');
c_data_set_4x = getappdata(0,'c_dataset_4x');
c_data_set_4y = getappdata(0,'c_dataset_4y');
c_data_set_4d = getappdata(0,'c_dataset_4d');


l_spec_1 = getappdata(0,'l_spec1');
l_spec_2 = getappdata(0,'l_spec2');
l_spec_3 = getappdata(0,'l_spec3');
l_spec_4 = getappdata(0,'l_spec4');
length1 = getappdata(0,'length1');

%get tube radius data
s1_r = getappdata(0,'s1r');
s2_r = getappdata(0,'s2r');
s3_r = getappdata(0,'s3r');
s4_r = getappdata(0,'s4r');



fs = 34000*11/(2*length1);
x = (1:512)*fs/(2*512);

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        [pks,loc]=findpeaks(log(l_spec_1));
        F1 = x(loc(1));
        F2 = x(loc(2));
        F3 = x(loc(3));
        set(handles.F1, 'String', F1);
        set(handles.F2, 'String', F2);
        set(handles.F3, 'String', F3);
        plot(handles.axes1, data_set_1x, data_set_1y, 'red', data_set_1x, (data_set_1y+data_set_1d),data_set_1x, (data_set_1y-data_set_1d), 'blue' );
        plot(handles.axes5,c_data_set_1x, c_data_set_1y, 'red', c_data_set_4x, (c_data_set_1y+c_data_set_1d),c_data_set_1x, (c_data_set_1y-c_data_set_1d), 'blue');
        axes(handles.axes8);
        bar(s1_r,1);
        hold on
        bar(-1*s1_r,1);
        hold off
        xlim([0 12])
        axes(handles.axes9);
        plot (x, log(l_spec_1))
        set(handles.slider3, 'Value', s1_r(1));
        set(handles.slider_1, 'String', s1_r(1));
        set(handles.slider4, 'Value', s1_r(2));
        set(handles.slider_2, 'String', s1_r(2));
        set(handles.slider5, 'Value', s1_r(3));
        set(handles.slider_3, 'String', s1_r(3));
        set(handles.slider6, 'Value', s1_r(4));
        set(handles.slider_4, 'String', s1_r(4));
        set(handles.slider7, 'Value', s1_r(5));
        set(handles.slider_5, 'String', s1_r(5));
        set(handles.slider8, 'Value', s1_r(6));
        set(handles.slider_6, 'String', s1_r(6));
        set(handles.slider9, 'Value', s1_r(7));
        set(handles.slider_7, 'String', s1_r(7));
        set(handles.slider10, 'Value', s1_r(8));
        set(handles.slider_8, 'String', s1_r(8));
        set(handles.slider11, 'Value', s1_r(9));
        set(handles.slider_9, 'String', s1_r(9));
        set(handles.slider12, 'Value', s1_r(10));
        set(handles.slider_10, 'String', s1_r(10));
        set(handles.slider13, 'Value', s1_r(11));
        set(handles.slider_11, 'String', s1_r(11));
    case 2
        [pks,loc]=findpeaks(log(l_spec_2));
        F1 = x(loc(1));
        F2 = x(loc(2));
        F3 = x(loc(3));
        set(handles.F1, 'String', F1);
        set(handles.F2, 'String', F2);
        set(handles.F3, 'String', F3);
        plot(handles.axes1, data_set_2x, data_set_2y, 'red', data_set_2x, (data_set_2y+data_set_2d),data_set_2x, (data_set_2y-data_set_2d), 'blue' );
        plot(handles.axes5, c_data_set_2x, c_data_set_2y, 'red', c_data_set_4x, (c_data_set_2y+c_data_set_2d),c_data_set_2x, (c_data_set_2y-c_data_set_2d), 'blue');
        axes(handles.axes8);
        bar(s2_r,1);
        hold on
        bar(-1*s2_r,1);
        hold off
        xlim([0 12])
        axes(handles.axes9);
        plot (x, log(l_spec_2))
        
        set(handles.slider3, 'Value', s2_r(1));
        set(handles.slider_1, 'String', s2_r(1));
        set(handles.slider4, 'Value', s2_r(2));
        set(handles.slider_2, 'String', s2_r(2));
        set(handles.slider5, 'Value', s2_r(3));
        set(handles.slider_3, 'String', s2_r(3));
        set(handles.slider6, 'Value', s2_r(4));
        set(handles.slider_4, 'String', s2_r(4));
        set(handles.slider7, 'Value', s2_r(5));
        set(handles.slider_5, 'String', s2_r(5));
        set(handles.slider8, 'Value', s2_r(6));
        set(handles.slider_6, 'String', s2_r(6));
        set(handles.slider9, 'Value', s2_r(7));
        set(handles.slider_7, 'String', s2_r(7));
        set(handles.slider10, 'Value', s2_r(8));
        set(handles.slider_8, 'String', s2_r(8));
        set(handles.slider11, 'Value', s2_r(9));
        set(handles.slider_9, 'String', s2_r(9));
        set(handles.slider12, 'Value', s2_r(10));
        set(handles.slider_10, 'String', s2_r(10));
        set(handles.slider13, 'Value', s2_r(11));
        set(handles.slider_11, 'String', s2_r(11));
    case 3
        [ ~,loc]=findpeaks(log(l_spec_3));
        F1 = x(loc(1));
        F2 = x(loc(2));
        F3 = x(loc(3));
        set(handles.F1, 'String', F1);
        set(handles.F2, 'String', F2);
        set(handles.F3, 'String', F3);
        plot(handles.axes1, data_set_3x, data_set_3y, 'red', data_set_3x, (data_set_3y+data_set_3d),data_set_3x, (data_set_3y-data_set_3d), 'blue' );
        plot(handles.axes5, c_data_set_3x, c_data_set_3y, 'red', c_data_set_3x, (c_data_set_3y+c_data_set_3d),c_data_set_3x, (c_data_set_3y-c_data_set_3d), 'blue');
        axes(handles.axes8);
        bar(s3_r,1);
        hold on
        bar(-1*s3_r,1);
        hold off
        xlim([0 12])
        axes(handles.axes9);
        plot (x, log(l_spec_3))
        
        set(handles.slider3, 'Value', s3_r(1));
        set(handles.slider_1, 'String', s3_r(1));
        set(handles.slider4, 'Value', s3_r(2));
        set(handles.slider_2, 'String', s3_r(2));
        set(handles.slider5, 'Value', s3_r(3));
        set(handles.slider_3, 'String', s3_r(3));
        set(handles.slider6, 'Value', s3_r(4));
        set(handles.slider_4, 'String', s3_r(4));
        set(handles.slider7, 'Value', s3_r(5));
        set(handles.slider_5, 'String', s3_r(5));
        set(handles.slider8, 'Value', s3_r(6));
        set(handles.slider_6, 'String', s3_r(6));
        set(handles.slider9, 'Value', s3_r(7));
        set(handles.slider_7, 'String', s3_r(7));
        set(handles.slider10, 'Value', s3_r(8));
        set(handles.slider_8, 'String', s3_r(8));
        set(handles.slider11, 'Value', s3_r(9));
        set(handles.slider_9, 'String', s3_r(9));
        set(handles.slider12, 'Value', s3_r(10));
        set(handles.slider_10, 'String', s3_r(10));
        set(handles.slider13, 'Value', s3_r(11));
        set(handles.slider_11, 'String', s3_r(11));
    case 4
        [pks,loc]=findpeaks(log(l_spec_4));
        F1 = x(loc(1));
        F2 = x(loc(2));
        F3 = x(loc(3));
        set(handles.F1, 'String', F1);
        set(handles.F2, 'String', F2);
        set(handles.F3, 'String', F3);
        plot(handles.axes1, data_set_4x, data_set_4y, 'red', data_set_4x, (data_set_4y+data_set_4d),data_set_4x, (data_set_4y-data_set_4d), 'blue' );
        plot(handles.axes5, c_data_set_4x, c_data_set_4y, 'red', c_data_set_4x, (c_data_set_4y+c_data_set_4d),c_data_set_4x, (c_data_set_4y-c_data_set_4d), 'blue');
        axes(handles.axes8);
        bar(s4_r,1);
        hold on
        bar(-1*s4_r,1);
        hold off
        xlim([0 12])
        axes(handles.axes9);
        plot (x, log(l_spec_4))
        
        set(handles.slider3, 'Value', s4_r(1));
        set(handles.slider_1, 'String', s4_r(1));
        set(handles.slider4, 'Value', s4_r(2));
        set(handles.slider_2, 'String', s4_r(2));
        set(handles.slider5, 'Value', s4_r(3));
        set(handles.slider_3, 'String', s4_r(3));
        set(handles.slider6, 'Value', s4_r(4));
        set(handles.slider_4, 'String', s4_r(4));
        set(handles.slider7, 'Value', s4_r(5));
        set(handles.slider_5, 'String', s4_r(5));
        set(handles.slider8, 'Value', s4_r(6));
        set(handles.slider_6, 'String', s4_r(6));
        set(handles.slider9, 'Value', s4_r(7));
        set(handles.slider_7, 'String', s4_r(7));
        set(handles.slider10, 'Value', s4_r(8));
        set(handles.slider_8, 'String', s4_r(8));
        set(handles.slider11, 'Value', s4_r(9));
        set(handles.slider_9, 'String', s4_r(9));
        set(handles.slider12, 'Value', s4_r(10));
        set(handles.slider_10, 'String', s4_r(10));
        set(handles.slider13, 'Value', s4_r(11));
        set(handles.slider_11, 'String', s4_r(11));
end






% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
    ['Close ' get(handles.figure1,'Name') '...'],...
    'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
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

set(hObject, 'String', {'Set 1', 'Set 2', 'Set 3', 'Set 4'});


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

start = 13;


set(handles.start, 'String', start);

popup_sel_index = get(handles.popupmenu2, 'Value');
switch popup_sel_index
    case 1
        datatype = 0;
    case 2
        datatype = 1;
        
end

Data_all = sel_file(datatype);
setappdata(0,'Data_all',Data_all);

if datatype == 0
Data = Data_all(1);
else
Data = Data_all{1};
end
tic;
%data as obtained from the results of the scanning
data_set_1x = ext_data(Data, datatype, 1);
data_set_1y = ext_data(Data, datatype, 2);
data_set_1d = ext_data(Data, datatype, 3);
data_set_2x = ext_data(Data, datatype, 4);
data_set_2y = ext_data(Data, datatype, 5);
data_set_2d = ext_data(Data, datatype, 6);
data_set_3x = ext_data(Data, datatype, 7);
data_set_3y = ext_data(Data, datatype, 8);
data_set_3d = ext_data(Data, datatype, 9);
data_set_4x = ext_data(Data, datatype, 10);
data_set_4y = ext_data(Data, datatype, 11);
data_set_4d = ext_data(Data, datatype, 12);



%find glottis

if datatype == 0
    finish = find_finish(data_set_1y);
else
    finish = 55;
end
setappdata(0,'finish',finish);
set(handles.finish, 'String', finish);

%excess data from area outside the region of interest removed

c_data_set_1x = ext_ext_data(data_set_1x, start, finish, datatype);
c_data_set_1y = ext_ext_data(data_set_1y, start, finish, datatype);
c_data_set_1d = ext_ext_data(data_set_1d, start, finish, datatype);
c_data_set_2x = ext_ext_data(data_set_2x, start, finish, datatype);
c_data_set_2y = ext_ext_data(data_set_2y, start, finish, datatype);
c_data_set_2d = ext_ext_data(data_set_2d, start, finish, datatype);
c_data_set_3x = ext_ext_data(data_set_3x, start, finish, datatype);
c_data_set_3y = ext_ext_data(data_set_3y, start, finish, datatype);
c_data_set_3d = ext_ext_data(data_set_3d, start, finish, datatype);
c_data_set_4x = ext_ext_data(data_set_4x, start, finish, datatype);
c_data_set_4y = ext_ext_data(data_set_4y, start, finish, datatype);
c_data_set_4d = ext_ext_data(data_set_4d, start, finish, datatype);

%set data for other functions to use
setappdata(0,'dataset_1x',data_set_1x);
setappdata(0,'dataset_1y',data_set_1y);
setappdata(0,'dataset_1d',data_set_1d);
setappdata(0,'dataset_2x',data_set_2x);
setappdata(0,'dataset_2y',data_set_2y);
setappdata(0,'dataset_2d',data_set_2d);
setappdata(0,'dataset_3x',data_set_3x);
setappdata(0,'dataset_3y',data_set_3y);
setappdata(0,'dataset_3d',data_set_3d);
setappdata(0,'dataset_4x',data_set_4x);
setappdata(0,'dataset_4y',data_set_4y);
setappdata(0,'dataset_4d',data_set_4d);

setappdata(0,'c_dataset_1x',c_data_set_1x);
setappdata(0,'c_dataset_1y',c_data_set_1y);
setappdata(0,'c_dataset_1d',c_data_set_1d);
setappdata(0,'c_dataset_2x',c_data_set_2x);
setappdata(0,'c_dataset_2y',c_data_set_2y);
setappdata(0,'c_dataset_2d',c_data_set_2d);
setappdata(0,'c_dataset_3x',c_data_set_3x);
setappdata(0,'c_dataset_3y',c_data_set_3y);
setappdata(0,'c_dataset_3d',c_data_set_3d);
setappdata(0,'c_dataset_4x',c_data_set_4x);
setappdata(0,'c_dataset_4y',c_data_set_4y);
setappdata(0,'c_dataset_4d',c_data_set_4d);


if datatype == 0
    
    %set up lenth
    length1 = data_set_1x(finish) - data_set_1x(start);
    length2 = data_set_2x(finish) - data_set_2x(start);
    length3 = data_set_3x(finish) - data_set_3x(start);
    length4 = data_set_4x(finish) - data_set_4x(start);
    %set up tube length
    tubelength2 = length2/11;
    tubelength3 = length3/11;
    tubelength4 = length4/11;
    %radius
    s2_r = get_tubes_average(data_set_2x, data_set_2y, tubelength2, start, finish, datatype);
    s3_r = get_tubes_average(data_set_3x, data_set_3y, tubelength3, start, finish, datatype);
    s4_r = get_tubes_average(data_set_4x, data_set_4y, tubelength4, start, finish, datatype);
    %spectrum
    l_spec_2 = f_spec(s2_r);
    l_spec_3 = f_spec(s3_r);
    l_spec_4 = f_spec(s4_r);
    %set data
    setappdata(0,'s2r',s2_r);
    setappdata(0,'s3r',s3_r);
    setappdata(0,'s4r',s4_r);
    setappdata(0,'l_spec2',l_spec_2);
    setappdata(0,'l_spec3',l_spec_3);
    setappdata(0,'l_spec4',l_spec_4);
    
else
    length1 =  data_set_1x(size(data_set_1x,1));
end

tubelength1 = length1/11;
s1_r = get_tubes_average(data_set_1x, data_set_1y, tubelength1, start, finish, datatype);
l_spec_1 = f_spec(s1_r);
setappdata(0,'length1',length1);
setappdata(0,'s1r',s1_r);
setappdata(0,'l_spec1',l_spec_1);

toc;



function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start as text
%        str2double(get(hObject,'String')) returns contents of start as a double


% --- Executes during object creation, after setting all properties.
function start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function finish_Callback(hObject, eventdata, handles)
% hObject    handle to finish (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of finish as text
%        str2double(get(hObject,'String')) returns contents of finish as a double


% --- Executes during object creation, after setting all properties.
function finish_CreateFcn(hObject, eventdata, handles)
% hObject    handle to finish (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2

contents = cellstr(get(hObject,'String'));
selection = contents{get(hObject,'Value')};
if strcmp(selection, 'AR')==1
        set(handles.vowel_choice, 'String', {'EXHALE', 'GLOT', 'HAD', 'HARD', 'HEAD', 'HEED', 'HERD', 'HID', 'HOARD', 'HOD', 'WHOD'});
else
        set(handles.vowel_choice, 'String', {'HAD', 'HARD', 'HEAD', 'HEED', 'HERD', 'HID', 'HOARD', 'HOD', 'HOOD', 'HOD'});

        
end

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


set(hObject, 'String', {'AR', 'MRI'});


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

slider_value = get(handles.slider3,'Value');
slider_1 = slider_value;
set(handles.slider_1, 'String', slider_1);
set(handles.slider3, 'Value', slider_1);

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

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_value = get(handles.slider4,'Value');
slider_2 = slider_value;
set(handles.slider_2, 'String', slider_2);
set(handles.slider4, 'Value', slider_2);

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

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_value = get(handles.slider5,'Value');
slider_3 = slider_value;
set(handles.slider_3, 'String', slider_3);
set(handles.slider5, 'Value', slider_3);

% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_value = get(handles.slider6,'Value');
slider_4 = slider_value;
set(handles.slider_4, 'String', slider_4);
set(handles.slider6, 'Value', slider_4);

% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider7_Callback(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_value = get(handles.slider7,'Value');
slider_5 = slider_value;
set(handles.slider_5, 'String', slider_5);
set(handles.slider7, 'Value', slider_5);

% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider8_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_value = get(handles.slider8,'Value');
slider_6 = slider_value;
set(handles.slider_6, 'String', slider_6);
set(handles.slider8, 'Value', slider_6);

% --- Executes during object creation, after setting all properties.
function slider8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider9_Callback(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_value = get(handles.slider9,'Value');
slider_7 = slider_value;
set(handles.slider_7, 'String', slider_7);
set(handles.slider9, 'Value', slider_7);

% --- Executes during object creation, after setting all properties.
function slider9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider10_Callback(hObject, eventdata, handles)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_value = get(handles.slider10,'Value');
slider_8 = slider_value;
set(handles.slider_8, 'String', slider_8);
set(handles.slider10, 'Value', slider_8);

% --- Executes during object creation, after setting all properties.
function slider10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider11_Callback(hObject, eventdata, handles)
% hObject    handle to slider11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of
%        slider
slider_value = get(handles.slider11,'Value');
slider_9 = slider_value;
set(handles.slider_9, 'String', slider_9);
set(handles.slider11, 'Value', slider_9);

% --- Executes during object creation, after setting all properties.
function slider11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider12_Callback(hObject, eventdata, handles)
% hObject    handle to slider12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_value = get(handles.slider12,'Value');
slider_10 = slider_value;
set(handles.slider_10, 'String', slider_10);
set(handles.slider12, 'Value', slider_10);

% --- Executes during object creation, after setting all properties.
function slider12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider13_Callback(hObject, eventdata, handles)
% hObject    handle to slider13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
slider_value = get(handles.slider13,'Value');
slider_11 = slider_value;
set(handles.slider_11, 'String', slider_11);
set(handles.slider13, 'Value', slider_11);

% --- Executes during object creation, after setting all properties.
function slider13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider17_Callback(hObject, eventdata, handles)
% hObject    handle to slider17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

slider_value = get(handles.finish,'Value');
slider_12 = slider_value;
x(1:11) = slider_value;
y = 0:10;
plot(handles.axes5,x,y)
slider_value2 = get(handles.slider17,'Value');
set(handles.finish, 'Value', slider_value2);

% --- Executes during object creation, after setting all properties.
function slider17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in vowel_choice.
function vowel_choice_Callback(hObject, eventdata, handles)
% hObject    handle to vowel_choice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns vowel_choice contents as cell array
%        contents{get(hObject,'Value')} returns selected item from vowel_choice


start = 13;


set(handles.start, 'String', start);


Data_all = getappdata(0,'Data_all');

popup_sel_index = get(handles.popupmenu2, 'Value');
switch popup_sel_index
    case 1
        datatype = 0;

popup_sel_index = get(handles.vowel_choice, 'Value');
switch popup_sel_index
    case 1
Data = Data_all(1);        
    case 2
Data = Data_all(2);        
    case 3
Data = Data_all(3);        
    case 4
Data = Data_all(4);        
    case 5
Data = Data_all(5);        
    case 6
Data = Data_all(6);        
    case 7
Data = Data_all(7);        
    case 8
Data = Data_all(8);        
    case 9
Data = Data_all(9);        
    case 10
Data = Data_all(10);        
    case 11
Data = Data_all(11);        
end
    case 2
        datatype = 1;

popup_sel_index = get(handles.vowel_choice, 'Value');
switch popup_sel_index
    case 1
Data = Data_all{1};        
    case 2
Data = Data_all{2};        
    case 3
Data = Data_all{3};        
    case 4
Data = Data_all{4};        
    case 5
Data = Data_all{5};        
    case 6
Data = Data_all{6};        
    case 7
Data = Data_all{7};        
    case 8
Data = Data_all{8};        
    case 9
Data = Data_all{9};        
    case 10
Data = Data_all{10};            
end
end



%data as obtained from the results of the scanning
data_set_1x = ext_data(Data, datatype, 1);
data_set_1y = ext_data(Data, datatype, 2);
data_set_1d = ext_data(Data, datatype, 3);
data_set_2x = ext_data(Data, datatype, 4);
data_set_2y = ext_data(Data, datatype, 5);
data_set_2d = ext_data(Data, datatype, 6);
data_set_3x = ext_data(Data, datatype, 7);
data_set_3y = ext_data(Data, datatype, 8);
data_set_3d = ext_data(Data, datatype, 9);
data_set_4x = ext_data(Data, datatype, 10);
data_set_4y = ext_data(Data, datatype, 11);
data_set_4d = ext_data(Data, datatype, 12);



%find glottis

if datatype == 0
    finish = find_finish(data_set_1y);
else
    finish = 55;
end
setappdata(0,'finish',finish);
set(handles.finish, 'String', finish);

%excess data from area outside the region of interest removed

c_data_set_1x = ext_ext_data(data_set_1x, start, finish, datatype);
c_data_set_1y = ext_ext_data(data_set_1y, start, finish, datatype);
c_data_set_1d = ext_ext_data(data_set_1d, start, finish, datatype);
c_data_set_2x = ext_ext_data(data_set_2x, start, finish, datatype);
c_data_set_2y = ext_ext_data(data_set_2y, start, finish, datatype);
c_data_set_2d = ext_ext_data(data_set_2d, start, finish, datatype);
c_data_set_3x = ext_ext_data(data_set_3x, start, finish, datatype);
c_data_set_3y = ext_ext_data(data_set_3y, start, finish, datatype);
c_data_set_3d = ext_ext_data(data_set_3d, start, finish, datatype);
c_data_set_4x = ext_ext_data(data_set_4x, start, finish, datatype);
c_data_set_4y = ext_ext_data(data_set_4y, start, finish, datatype);
c_data_set_4d = ext_ext_data(data_set_4d, start, finish, datatype);

%set data for other functions to use
setappdata(0,'dataset_1x',data_set_1x);
setappdata(0,'dataset_1y',data_set_1y);
setappdata(0,'dataset_1d',data_set_1d);
setappdata(0,'dataset_2x',data_set_2x);
setappdata(0,'dataset_2y',data_set_2y);
setappdata(0,'dataset_2d',data_set_2d);
setappdata(0,'dataset_3x',data_set_3x);
setappdata(0,'dataset_3y',data_set_3y);
setappdata(0,'dataset_3d',data_set_3d);
setappdata(0,'dataset_4x',data_set_4x);
setappdata(0,'dataset_4y',data_set_4y);
setappdata(0,'dataset_4d',data_set_4d);

setappdata(0,'c_dataset_1x',c_data_set_1x);
setappdata(0,'c_dataset_1y',c_data_set_1y);
setappdata(0,'c_dataset_1d',c_data_set_1d);
setappdata(0,'c_dataset_2x',c_data_set_2x);
setappdata(0,'c_dataset_2y',c_data_set_2y);
setappdata(0,'c_dataset_2d',c_data_set_2d);
setappdata(0,'c_dataset_3x',c_data_set_3x);
setappdata(0,'c_dataset_3y',c_data_set_3y);
setappdata(0,'c_dataset_3d',c_data_set_3d);
setappdata(0,'c_dataset_4x',c_data_set_4x);
setappdata(0,'c_dataset_4y',c_data_set_4y);
setappdata(0,'c_dataset_4d',c_data_set_4d);


if datatype == 0
    
    %set up lenth
    length1 = data_set_1x(finish) - data_set_1x(start);
    length2 = data_set_2x(finish) - data_set_2x(start);
    length3 = data_set_3x(finish) - data_set_3x(start);
    length4 = data_set_4x(finish) - data_set_4x(start);
    %set up tube length
    tubelength2 = length2/11;
    tubelength3 = length3/11;
    tubelength4 = length4/11;
    %radius
    s2_r = get_tubes_average(data_set_2x, data_set_2y, tubelength2, start, finish, datatype);
    s3_r = get_tubes_average(data_set_3x, data_set_3y, tubelength3, start, finish, datatype);
    s4_r = get_tubes_average(data_set_4x, data_set_4y, tubelength4, start, finish, datatype);
    %spectrum
    l_spec_2 = f_spec(s2_r);
    l_spec_3 = f_spec(s3_r);
    l_spec_4 = f_spec(s4_r);
    %set data
    setappdata(0,'s2r',s2_r);
    setappdata(0,'s3r',s3_r);
    setappdata(0,'s4r',s4_r);
    setappdata(0,'l_spec2',l_spec_2);
    setappdata(0,'l_spec3',l_spec_3);
    setappdata(0,'l_spec4',l_spec_4);
    
else
    length1 =  data_set_1x(size(data_set_1x,1));
end

tubelength1 = length1/11;
s1_r = get_tubes_average(data_set_1x, data_set_1y, tubelength1, start, finish, datatype);
l_spec_1 = f_spec(s1_r);
setappdata(0,'length1',length1);
setappdata(0,'s1r',s1_r);
setappdata(0,'l_spec1',l_spec_1);

% --- Executes during object creation, after setting all properties.
function vowel_choice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vowel_choice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
