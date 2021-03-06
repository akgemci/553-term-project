%% AYSEL AKGEMC� - 1875566 
% Fall 2016 Term Project
% OPTIMIZATION OF THE MASS OF A RFPM GENERATOR USING PSO
%% ----------------------- MAIN--------------------------
% First, run the  main code which is AA_termproject.m 
% then push 'RUN' button on the GUI
function varargout = AA_termproject(varargin)
% AA_TERMPROJECT MATLAB code for AA_termproject.fig
%      AA_TERMPROJECT, by itself, creates a new AA_TERMPROJECT or raises the existing
%      singleton*.
%
%      H = AA_TERMPROJECT returns the handle to a new AA_TERMPROJECT or the handle to
%      the existing singleton*.
%
%      AA_TERMPROJECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AA_TERMPROJECT.M with the given input arguments.
%
%      AA_TERMPROJECT('Property','Value',...) creates a new AA_TERMPROJECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AA_termproject_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AA_termproject_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AA_termproject

% Last Modified by GUIDE v2.5 29-Jan-2017 15:57:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AA_termproject_OpeningFcn, ...
                   'gui_OutputFcn',  @AA_termproject_OutputFcn, ...
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


% --- Executes just before AA_termproject is made visible.
function AA_termproject_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AA_termproject (see VARARGIN)

% Choose default command line output for AA_termproject
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AA_termproject wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AA_termproject_OutputFcn(hObject, eventdata, handles) 
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

x = str2double(get(handles.iternum,'string'));
results = zeros(1,4);
[results(1) results(2) results(3) results(4)] = PSO(x);
set(handles.lm, 'string', round(results(1)*10000)/10);
set(handles.di, 'string', round(results(2)*10000)/10);
set(handles.pp, 'string', round(results(3)));
set(handles.mass, 'string', results(4));



function iternum_Callback(hObject, eventdata, handles)
% hObject    handle to iternum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iternum as text
%        str2double(get(hObject,'String')) returns contents of iternum as a double


% --- Executes during object creation, after setting all properties.
function iternum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iternum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
