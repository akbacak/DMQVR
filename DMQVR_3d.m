function varargout = DMQVR_3d(varargin)
% DMQVR_3D MATLAB code for DMQVR_3d.fig
%      DMQVR_3D, by itself, creates a new DMQVR_3D or raises the existing
%      singleton*.
%
%      H = DMQVR_3D returns the handle to a new DMQVR_3D or the handle to
%      the existing singleton*.
%
%      DMQVR_3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DMQVR_3D.M with the given input arguments.
%
%      DMQVR_3D('Property','Value',...) creates a new DMQVR_3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DMQVR_3d_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DMQVR_3d_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DMQVR_3d

% Last Modified by GUIDE v2.5 06-May-2019 23:44:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DMQVR_3d_OpeningFcn, ...
                   'gui_OutputFcn',  @DMQVR_3d_OutputFcn, ...
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


% --- Executes just before DMQVR_3d is made visible.
function DMQVR_3d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DMQVR_3d (see VARARGIN)

% Choose default command line output for DMQVR_3d
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DMQVR_3d wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DMQVR_3d_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Dataset.
function Dataset_Callback(hObject, eventdata, handles)
% hObject    handle to Dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Dataset contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Dataset

maxFront = 3;

dataset_index = get(handles.Dataset, 'Value');
switch dataset_index 
    case 1
        video_dir =[pwd '/myDataset/videoFolder/']; 
        data_dir = [pwd '/myDataset/hashCodes/'];
        
        colorData = 0;   
        
    case 2
        video_dir =[pwd '/myDataset2/videoFolder/']; 
        data_dir = [pwd '/myDataset2/hashCodes/'];
       
        colorData = 0;    
end
% if(colorData == 1)
% load([data_dir '/colorLabel']);
% handles.colorLabel = colorLabel;
% end


load([data_dir '/filenames']);
load([data_dir '/targets']);


handles.filenames = filenames;
handles.targets   = targets;
handles.video_dir = video_dir;
handles.data_dir  = data_dir;
handles.maxFront  = maxFront;


set(handles.QueryName1,'String', filenames);
set(handles.QueryName2,'String', filenames);
set(handles.QueryName3,'String', filenames);



guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Dataset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dataset_index = get(handles.Dataset, 'Value');
switch dataset_index 
    case 1
        image_dir =[pwd '/lamdaDataset/scene_categories/']; 
        data_dir = [pwd '/lamdaDataset/preprocessed_features'];
        colorData = 0;   
        
    case 2
        image_dir =[pwd '/newdbDataset/scene_categories/']; 
        data_dir = [pwd '/newdbDataset/preprocessed_features'];
        colorData = 0;    
end
% if(colorData == 1)
% load([data_dir '/colorLabel']);
% handles.colorLabel = colorLabel;
% end

set(handles.Status,'String','Loading...');pause(0.3);
load([data_dir '/filenames']);
load([data_dir '/targets']);
load([data_dir '/pyramid_all']);set(handles.Status,'String','Loading 50%...');pause(0.3);

handles.filenames = filenames;
handles.data = pyramid_all;
handles.targets = targets;

set(handles.Status,'String','Done');
set(handles.QueryName1,'String', filenames);
set(handles.QueryName2,'String', filenames);
set(handles.QueryName3,'String', filenames);

handles.image_dir = image_dir;
set(handles.Status,'String','Loading');
set(handles.Status,'String','Done');

guidata(hObject, handles);


% --- Executes on selection change in QueryName1.
function QueryName1_Callback(hObject, eventdata, handles)
% hObject    handle to QueryName1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns QueryName1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from QueryName1

imindex = get(hObject,'Value');
video_dir = handles.video_dir;
fname = [video_dir handles.filenames{imindex}];

mov=VideoReader(fname);
nFrames=mov.NumberOfFrames;
for i=1:nFrames
  videoFrame=read(mov,i);
  axes(handles.axes1);
  imshow(videoFrame);  axis image;
end




handles.q1Idx = imindex;
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function QueryName1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to QueryName1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in QueryName2.
function QueryName2_Callback(hObject, eventdata, handles)
% hObject    handle to QueryName2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns QueryName2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from QueryName2
imindex = get(hObject,'Value');
video_dir = handles.video_dir;
fname = [video_dir handles.filenames{imindex}];

mov=VideoReader(fname);
nFrames=mov.NumberOfFrames;
for i=1:nFrames
  videoFrame=read(mov,i);
  axes(handles.axes2);
  imshow(videoFrame);  axis image;
end


handles.q2Idx = imindex;
guidata(hObject, handles);


% --- Executes on selection change in QueryName3.
function QueryName3_Callback(hObject, eventdata, handles)
% hObject    handle to QueryName3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns QueryName3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from QueryName3
imindex = get(hObject,'Value');
video_dir = handles.video_dir;
fname = [video_dir handles.filenames{imindex}];

mov=VideoReader(fname);
nFrames=mov.NumberOfFrames;
for i=1:nFrames
  videoFrame=read(mov,i);
  axes(handles.axes12);
  imshow(videoFrame);  axis image;
end


handles.q3Idx = imindex;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function QueryName2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to QueryName2 (see GCBO)
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
    
    targets  = handles.targets;
    maxFront = handles.maxFront;

    queryIndex1 = handles.q1Idx;
    queryIndex2 = handles.q2Idx;
    queryIndex3 = handles.q3Idx;
    
    data = handles.data;
    q1 = data(queryIndex1,:);
    q2 = data(queryIndex2,:);
    q3 = data(queryIndex3,:);
    
      
    N = length(handles.filenames);  % N = 2000
 tic
    q1new = repmat(q1,N,1);
    q2new = repmat(q2,N,1);
    q3new = repmat(q3,N,1);
    dist_1 = xor(data, q1new);
    dist_2 = xor(data, q2new);
    dist_3 = xor(data, q3new);
    hamming_dist1 = sum(dist_1,2);
    hamming_dist2 = sum(dist_2,2);
    hamming_dist3 = sum(dist_3,2);
    n_hamming_dist1 = mat2gray(hamming_dist1);
    n_hamming_dist2 = mat2gray(hamming_dist2);
    n_hamming_dist3 = mat2gray(hamming_dist3);
           
t = toc;
set(handles.tictoc,'String',num2str(t))

    X = zeros(3,N);
    X(1,:) = n_hamming_dist1;
    X(2,:) = n_hamming_dist2;
    X(3,:) = n_hamming_dist3;
    
    X = (X)';
    [K,L] = size(unique(X,'rows'));  %% Number of unique pareto points 
    set(handles.num_pp_3d,'String',num2str(K))

    axes(handles.axes3);
    hold off; 
    scatter3(X(:,1),X(:,2),X(:,3),'k.');
    %plot3(X(:,1),X(:,2),X(:,3),'.');
    hold on; 
    
     %scatter3(X(:,1),X(:,2),X(:,3),'k.');
     
    
   
    [pf_idx] = pareto_fronts(X, maxFront);
    for k=1:1
        scatter3(pf_idx{k,1}(:,1), pf_idx{k,1}(:,2) , pf_idx{k,1}(:,3)); 
        %plot(pf_idx{k,1}(:,1), pf_idx{k,1}(:,2) , 'y-');
    end
    xlabel('c1');
    ylabel('c2'); 

   
    handles.pf_idx = pf_idx;
    handles.X = X;
       
    
    guidata(hObject, handles);
    


% --- Executes on slider movement.
function FrontSelector_Callback(hObject, eventdata, handles)
% hObject    handle to FrontSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

targets  = handles.targets;
maxFront = handles.maxFront;


currentFront = ((round(1+(maxFront-1)*get(hObject,'Value'))));

%set(handles.FrontIdx,'String',['Pareto depth:' num2str(currentFront)] );
set(handles.FrontNum,'String',num2str(currentFront));



pf_idx = handles.pf_idx;
X = handles.X;

axes(handles.axes3);
hold off; scatter3(handles.X(:,1),handles.X(:,2),handles.X(:,3),'k.');
view(3);
rotate3d on;
hold on;

 
 [pf_idx] = pareto_fronts(X, maxFront);
 
 %for k=1:maxFront
 %       scatter3(pf_idx{k,1}(:,1), pf_idx{k,1}(:,2) , pf_idx{k,1}(:,3)); 
 %       %plot(pf_idx{k,1}(:,1), pf_idx{k,1}(:,2) , 'y-');
 %end

 
  l = currentFront; 
   scatter3(pf_idx{l,1}(:,1), pf_idx{l,1}(:,2) , pf_idx{l,1}(:,3), 'o', 'filled'); 
   view(3);
   rotate3d on;
  %plot(pf_idx{l,1}(:,1), pf_idx{l,1}(:,2) , 'b-*');
  xlabel('c1');
  ylabel('c2');
 

 handles.currentFront = currentFront;
    
    
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    queryIndex1 = handles.q1Idx;
    queryIndex2 = handles.q2Idx;
    queryIndex3 = handles.q3Idx;
    data = handles.data;
    q1 = data(queryIndex1,:);
    q2 = data(queryIndex2,:);
    q3 = data(queryIndex3,:);

    q1_label = targets(queryIndex1,: ); % Label vector of Query 1
    q2_label = targets(queryIndex2,: ); % Label vector of Query 2
    q3_label = targets(queryIndex3,: ); % Label vector of Query 2
        
    b = or(q1_label , q2_label); % beta in the equation 7  
    b = or(b,q3_label);
       
    absolute_b = nnz(b);         % Number of non-zero elements in the beta, nnz is a Matlab Func.
    
    %absolute_b = sum(q1_label,2) + sum(q2_label,2) + sum(q3_label,2); % beta in the equation 7
    
    
      for j  =1:maxFront
            R=0;
            C=0;
            Labels = targets(pf_idx{j,1}(:,4),:);     
            [R , C] = size(Labels); 
            
            
              switch (mod(R ,2) == 1)
                case 1
                     e_left(j) = (round(R / 2) - 1) ;  
                     e_rigth(j)= (round(R / 2)    ) ;
                case 0
                     e_left(j)  =  R / 2 ;
                     e_rigth(j) =  R / 2 ;
              end
            
              
               
                     
            for e = 1:R              
                MQUR_ALL(j,e) =  nnz( and(Labels(e,:) , b ) ) /  absolute_b;
                                           
            end                   
            %MQUR score of the first retrieved items in each front
            MQUR_1(j) = MQUR_ALL(j,1);
              
            G{j} = 0;
            rigth_G{j} = 0;
            left_G{j}  = 0;
            i_G{j} = 0;
            i_rigth_G{j} = 0;
            i_left_G{j}  = 0;
            
            rigth_DG{j} = 0;
            left_DG{j}  = 0;
            i_DG{j} = 0;
            i_rigth_DG{j} = 0;
            i_left_DG{j}  = 0;
            
            rigth_DCG{j} = 0;
            left_DCG{j}  = 0;
            i_rigth_DCG{j} = 0;
            i_left_DCG{j}  = 0;        
              
            
            
              for s = 1:R     
               G{j}(:,s) = MQUR_ALL(j,s);
               i_G{j}(:,s) = 1;
              end
             
              
            if R <= 2
                    e_left(j) = 2;
                    e_rigth(j)= 2;
            end
            
            
             switch (mod(R ,2) == 1)
                case 1                    
                    rigth_G{j}   =   G{j}(: , (e_rigth(j)   )  : R);
                    i_rigth_G{j} = i_G{j}(: , (e_rigth(j)   )  : R);
                    
                    left_G{j}    =   G{j}(: , 1 : e_left(j) ); 
                    left_G{j}    =   fliplr(left_G{j});
                    i_left_G{j}   = i_G{j}(: , 1 : e_left(j) );
             
                    rigth_DG{j}(:,1)   =   rigth_G{j}(:,1);
                    i_rigth_DG{j}(:,1) = i_rigth_G{j}(:,1);
                    
             
             
                    left_DG{j}(:,1)   =   left_G{j}(:,1);
                    i_left_DG{j}(:,1) = i_left_G{j}(:,1);
             
                    for d = 2 : e_rigth(j)
                         rigth_DG{j}(:,d)   =   rigth_G{j}(:,d) / log2(d);
                         i_rigth_DG{j}(:,d) = i_rigth_G{j}(:,d) / log2(d); 
                    end                    
                    for  f = 2 : e_left(j)   
                         left_DG{j}(:,f)   =   left_G{j}(:,f) / log2(f);
                         i_left_DG{j}(:,f) = i_left_G{j}(:,f) / log2(f);
                    end
                    
                    case 0
                    rigth_G{j}   =   G{j}(: , (e_rigth(j) +1)  : R);
                    i_rigth_G{j} = i_G{j}(: , (e_rigth(j) +1)  : R);
                    
                    left_G{j}    =   G{j}(: , 1 : e_left(j) ); 
                    left_G{j}    =   fliplr(left_G{j});
                    i_left_G{j}   = i_G{j}(: , 1 : e_left(j) );
             
                    rigth_DG{j}(:,1)   =   rigth_G{j}(:,1);
                    i_rigth_DG{j}(:,1) = i_rigth_G{j}(:,1);
                    %rigth_DG{l,j}(:,e_rigth(l,j)) =  MQUR_ALL{l,:}(j,e_rigth(l,j));
             
             
                    left_DG{j}(:,1)   =   left_G{j}(:,1);
                    i_left_DG{j}(:,1) = i_left_G{j}(:,1);
             
             
                    for h = 2 : e_rigth(j)
                         rigth_DG{j}(:,h)   =   rigth_G{j}(:,h) / log2(h);
                         i_rigth_DG{j}(:,h) = i_rigth_G{j}(:,h) / log2(h); 
                    end
                    for m = 2 : e_left(j)
                         left_DG{j}(:,m)   =   left_G{j}(:,m) / log2(m);
                         i_left_DG{j}(:,m) = i_left_G{j}(:,m) / log2(m);
                    end
             end
             
             rigth_DCG{j} = cumsum(rigth_DG{j}); 
             i_rigth_DCG{j} = cumsum( i_rigth_DG{j});
             n_rigth_DCG{j} = rigth_DCG{j} ./ i_rigth_DCG{j};
                                  
           
             left_DCG{j} = cumsum(left_DG{j}); 
             i_left_DCG{j} = cumsum( i_left_DG{j});
             n_left_DCG{j} = left_DCG{j} ./ i_left_DCG{j};
             
             flip_n_left_DCG{j} = fliplr(n_left_DCG{j});
             n_DCG{j} = horzcat(flip_n_left_DCG{j}, n_rigth_DCG{j});
             
              
            
      end
      
      j = currentFront;
      axes(handles.axes11);
      hold off;
      plot(n_DCG{j});
      ylim([0 1])
      xlabel('Number of Images');
      ylabel('nDCG');
      hold on;
      
      handles.currentFront = currentFront;
     
    
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    
    
    
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function FrontSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrontSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function ImageSelector_Callback(hObject, eventdata, handles)
% hObject    handle to ImageSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


maxFront = handles.maxFront;

pf_idx = handles.pf_idx;
X = handles.X;
currentFront = handles.currentFront ;

[frontSize, e] = size(pf_idx{currentFront,1}(:,4));
currentImage = ((round(1+(frontSize-1)*get(hObject,'Value'))));

set(handles.ImageNum,'String',num2str(currentImage));

%set(handles.PicIdx1,'String',num2str(max(currentImage-2,1)));
%set(handles.PicIdx2,'String',num2str(max(currentImage-1,1)));
%set(handles.PicIdx3,'String',num2str(currentImage));
%set(handles.PicIdx4,'String',num2str(min(currentImage+1,frontSize)));
%set(handles.PicIdx5,'String',num2str(min(currentImage+2,frontSize)));



L1 = max(1,currentImage - 1);
R1 = min(frontSize,currentImage+1);




axes(handles.axes3);
hold off; plot3(handles.X(:,1),handles.X(:,2),handles.X(:,3),'k.');
view(3)
rotate3d on;
hold on;



[pf_idx] = pareto_fronts(X, maxFront);
 %for k=1:maxFront
 %       plot3(pf_idx{k,1}(:,1), pf_idx{k,1}(:,2) ,pf_idx{k,1}(:,3), 'y-');
 %end


l = currentFront;

plot3(pf_idx{l,1}(:,1), pf_idx{l,1}(:,2) , pf_idx{l,1}(:,3) , 's');
view(3);
rotate3d on;

xlabel('Dissimilarity to Query1');
ylabel('Dissimilarity to Query2');
zlabel('Dissimilarity to Query3');

currentfrontLoc = (pf_idx{l,1}(:,1:3))';


plot3(currentfrontLoc(1,currentImage),currentfrontLoc(2,currentImage),currentfrontLoc(3,currentImage),'O','Linewidth',3,'MarkerSize',10);
plot3(currentfrontLoc(1,L1),currentfrontLoc(2,L1),currentfrontLoc(3,L1),'o','Linewidth',3,'MarkerSize',5);
plot3(currentfrontLoc(1,R1),currentfrontLoc(2,R1),currentfrontLoc(3,L1),'o','Linewidth',3,'MarkerSize',5);

    
currentfrontId = pf_idx{l,:}(:,4);
currentfrontId =(currentfrontId)';


  
%axes(handles.axes4);
fname = [handles.video_dir handles.filenames{currentfrontId(1,currentImage)}];
set(handles.main,'string',num2str(handles.filenames{currentfrontId(1,currentImage)}));
axis image;
mov=VideoReader(fname);
nFrames=mov.NumberOfFrames;
for i=1:10:nFrames
  videoFrame=read(mov,i);
  axes(handles.axes4);
  imshow(videoFrame);  axis image;
end

%axes(handles.axes5);
fname = [handles.video_dir handles.filenames{currentfrontId(1,R1)}];
set(handles.mainR1,'string',num2str(handles.filenames{currentfrontId(1,R1)}));
axis image;
mov=VideoReader(fname);
nFrames=mov.NumberOfFrames;
for i=1:nFrames
  videoFrame=read(mov,i);
  axes(handles.axes5);
  imshow(videoFrame);  axis image;
end

%axes(handles.axes8);
fname = [handles.video_dir handles.filenames{currentfrontId(1,L1)}];
set(handles.mainL1,'string',num2str(handles.filenames{currentfrontId(1,L1)}));
axis image;
mov=VideoReader(fname);
nFrames=mov.NumberOfFrames;
for i=1:nFrames
  videoFrame=read(mov,i);
  axes(handles.axes8);
  imshow(videoFrame);  axis image;
end



guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ImageSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function FrontNum_Callback(hObject, eventdata, handles)
% hObject    handle to FrontNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrontNum as text
%        str2double(get(hObject,'String')) returns contents of FrontNum as a double


% --- Executes during object creation, after setting all properties.
function FrontNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrontNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ImageNum_Callback(hObject, eventdata, handles)
% hObject    handle to ImageNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImageNum as text
%        str2double(get(hObject,'String')) returns contents of ImageNum as a double


% --- Executes during object creation, after setting all properties.
function ImageNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tictoc_Callback(hObject, eventdata, handles)
% hObject    handle to tictoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tictoc as text
%        str2double(get(hObject,'String')) returns contents of tictoc as a double


% --- Executes during object creation, after setting all properties.
function tictoc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tictoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function main_Callback(hObject, eventdata, handles)
% hObject    handle to main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of main as text
%        str2double(get(hObject,'String')) returns contents of main as a double


% --- Executes during object creation, after setting all properties.
function main_CreateFcn(hObject, eventdata, handles)
% hObject    handle to main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mainR1_Callback(hObject, eventdata, handles)
% hObject    handle to mainR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mainR1 as text
%        str2double(get(hObject,'String')) returns contents of mainR1 as a double


% --- Executes during object creation, after setting all properties.
function mainR1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainR1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mainR2_Callback(hObject, eventdata, handles)
% hObject    handle to mainR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mainR2 as text
%        str2double(get(hObject,'String')) returns contents of mainR2 as a double


% --- Executes during object creation, after setting all properties.
function mainR2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainR2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mainR3_Callback(hObject, eventdata, handles)
% hObject    handle to mainR3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mainR3 as text
%        str2double(get(hObject,'String')) returns contents of mainR3 as a double


% --- Executes during object creation, after setting all properties.
function mainR3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainR3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mainL1_Callback(hObject, eventdata, handles)
% hObject    handle to mainL1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mainL1 as text
%        str2double(get(hObject,'String')) returns contents of mainL1 as a double


% --- Executes during object creation, after setting all properties.
function mainL1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainL1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mainL2_Callback(hObject, eventdata, handles)
% hObject    handle to mainL2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mainL2 as text
%        str2double(get(hObject,'String')) returns contents of mainL2 as a double


% --- Executes during object creation, after setting all properties.
function mainL2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainL2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mainL3_Callback(hObject, eventdata, handles)
% hObject    handle to mainL3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mainL3 as text
%        str2double(get(hObject,'String')) returns contents of mainL3 as a double


% --- Executes during object creation, after setting all properties.
function mainL3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainL3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nDCG_Value_Callback(hObject, eventdata, handles)
% hObject    handle to nDCG_Value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nDCG_Value as text
%        str2double(get(hObject,'String')) returns contents of nDCG_Value as a double


% --- Executes during object creation, after setting all properties.
function nDCG_Value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nDCG_Value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataset_index = get(handles.Dataset, 'Value');
switch dataset_index 
    case 1
        image_dir =[pwd '/lamdaDataset/scene_categories/']; 
        data_dir = [pwd '/lamdaDataset/preprocessed_features'];
        colorData = 0;   
        
    case 2
        image_dir =[pwd '/newdbDataset/scene_categories/']; 
        data_dir = [pwd '/newdbDataset/preprocessed_features'];
        colorData = 0;    
end
% if(colorData == 1)
% load([data_dir '/colorLabel']);
% handles.colorLabel = colorLabel;
% end

set(handles.Status,'String','Loading...');pause(0.3);
load([data_dir '/filenames']);
load([data_dir '/targets']);
load([data_dir '/pyramid_all_48']);set(handles.Status,'String','Loading 50%...');pause(0.3);

handles.filenames = filenames;
handles.data = pyramid_all_48;
handles.targets = targets;

set(handles.Status,'String','Done');
set(handles.QueryName1,'String', filenames);
set(handles.QueryName2,'String', filenames);
set(handles.QueryName3,'String', filenames);

handles.image_dir = image_dir;
set(handles.Status,'String','Loading');
set(handles.Status,'String','Done');

guidata(hObject, handles);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataset_index = get(handles.Dataset, 'Value');
switch dataset_index 
    case 1
        image_dir =[pwd '/lamdaDataset/scene_categories/']; 
        data_dir = [pwd '/lamdaDataset/preprocessed_features'];
        colorData = 0;   
        
    case 2
        image_dir =[pwd '/newdbDataset/scene_categories/']; 
        data_dir = [pwd '/newdbDataset/preprocessed_features'];
        colorData = 0;    
end
% if(colorData == 1)
% load([data_dir '/colorLabel']);
% handles.colorLabel = colorLabel;
% end

set(handles.Status,'String','Loading...');pause(0.3);
load([data_dir '/filenames']);
load([data_dir '/targets']);
load([data_dir '/pyramid_all_64']);set(handles.Status,'String','Loading 50%...');pause(0.3);

handles.filenames = filenames;
handles.data = pyramid_all_64;
handles.targets = targets;

set(handles.Status,'String','Done');
set(handles.QueryName1,'String', filenames);
set(handles.QueryName2,'String', filenames);
set(handles.QueryName3,'String', filenames);

handles.image_dir = image_dir;
set(handles.Status,'String','Loading');
set(handles.Status,'String','Done');

guidata(hObject, handles);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataset_index = get(handles.Dataset, 'Value');
switch dataset_index 
    case 1
        image_dir =[pwd '/lamdaDataset/scene_categories/']; 
        data_dir = [pwd '/lamdaDataset/preprocessed_features'];
        colorData = 0;   
        
    case 2
        image_dir =[pwd '/newdbDataset/scene_categories/']; 
        data_dir = [pwd '/newdbDataset/preprocessed_features'];
        colorData = 0;    
end
% if(colorData == 1)
% load([data_dir '/colorLabel']);
% handles.colorLabel = colorLabel;
% end

set(handles.Status,'String','Loading...');pause(0.3);
load([data_dir '/filenames']);
load([data_dir '/targets']);
load([data_dir '/pyramid_all_96']);set(handles.Status,'String','Loading 50%...');pause(0.3);

handles.filenames = filenames;
handles.data = pyramid_all_96;
handles.targets = targets;

set(handles.Status,'String','Done');
set(handles.QueryName1,'String', filenames);
set(handles.QueryName2,'String', filenames);
set(handles.QueryName3,'String', filenames);

handles.image_dir = image_dir;
set(handles.Status,'String','Loading');
set(handles.Status,'String','Done');

guidata(hObject, handles);


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dataset_index = get(handles.Dataset, 'Value');
switch dataset_index 
    case 1
        image_dir =[pwd '/lamdaDataset/scene_categories/']; 
        data_dir = [pwd '/lamdaDataset/preprocessed_features'];
        colorData = 0;   
        
    case 2
        image_dir =[pwd '/newdbDataset/scene_categories/']; 
        data_dir = [pwd '/newdbDataset/preprocessed_features'];
        colorData = 0;    
end
% if(colorData == 1)
% load([data_dir '/colorLabel']);
% handles.colorLabel = colorLabel;
% end

set(handles.Status,'String','Loading...');pause(0.3);
load([data_dir '/filenames']);
load([data_dir '/targets']);
load([data_dir '/pyramid_all_256']);set(handles.Status,'String','Loading 50%...');pause(0.3);

handles.filenames = filenames;
handles.data = pyramid_all_256;
handles.targets = targets;

set(handles.Status,'String','Done');
set(handles.QueryName1,'String', filenames);
set(handles.QueryName2,'String', filenames);
set(handles.QueryName3,'String', filenames);

handles.image_dir = image_dir;
set(handles.Status,'String','Loading');
set(handles.Status,'String','Done');

guidata(hObject, handles);


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dataset_index = get(handles.Dataset, 'Value');
switch dataset_index 
    case 1
        image_dir =[pwd '/lamdaDataset/scene_categories/']; 
        data_dir = [pwd '/lamdaDataset/preprocessed_features'];
        colorData = 0;   
        
    case 2
        image_dir =[pwd '/newdbDataset/scene_categories/']; 
        data_dir = [pwd '/newdbDataset/preprocessed_features'];
        colorData = 0;    
end
% if(colorData == 1)
% load([data_dir '/colorLabel']);
% handles.colorLabel = colorLabel;
% end

set(handles.Status,'String','Loading...');pause(0.3);
load([data_dir '/filenames']);
load([data_dir '/targets']);
load([data_dir '/pyramid_all_Caffe_ssdh48_v3']);set(handles.Status,'String','Loading 50%...');pause(0.3);

handles.filenames = filenames;
handles.data = pyramid_all_Caffe_ssdh48_v3;
handles.targets = targets;

set(handles.Status,'String','Done');
set(handles.QueryName1,'String', filenames);
set(handles.QueryName2,'String', filenames);
set(handles.QueryName3,'String', filenames);

handles.image_dir = image_dir;
set(handles.Status,'String','Loading');
set(handles.Status,'String','Done');

guidata(hObject, handles);



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
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

 
    queryIndex1 = handles.q1Idx;
    queryIndex2 = handles.q2Idx;
    queryIndex3 = handles.q3Idx;
    data = handles.data;
    q1 = data(queryIndex1,:);
    q2 = data(queryIndex2,:);
    q3 = data(queryIndex3,:);
    
    
    %     Ranking with EMR
    N = length(handles.filenames);
tic    
    [H A landmarks Z] = EMRcomputeModel(handles.data);
    y1 = zeros(N,1);
    y1(queryIndex1) = 1;
    y2 = zeros(N,1);
    y2(queryIndex2) = 1;
    y3 = zeros(N,1);
    y3(queryIndex3) = 1;
    
    simEMR1 = EMRscore(H ,A, y1);
    simEMR2 = EMRscore(H ,A, y2);
    simEMR3 = EMRscore(H ,A, y3);
    dist1 = 1-simEMR1;
    dist2 = 1-simEMR2;
    dist3 = 1-simEMR3;
       
t = toc;
set(handles.tictoc2,'String',num2str(t))

    X = zeros(3,N);
    X(1,:) = dist1;
    X(2,:) = dist2;
    X(3,:) = dist3;
    
    X = (X)';
    
    [K,L] = size(unique(X,'rows'));  %% Number of unique pareto points 
    set(handles.num_pp_3d,'String',num2str(K))    
    
    axes(handles.axes3);
    hold off; 
    scatter3(X(:,1),X(:,2),X(:,3),'k.');
    %plot3(X(:,1),X(:,2),X(:,3),'.');
    hold on; 
    
     %scatter3(X(:,1),X(:,2),X(:,3),'k.');
     
    
  
    [pf_idx] = pareto_fronts(X, maxFront);
    for k=1:maxFront
        scatter3(pf_idx{k,1}(:,1), pf_idx{k,1}(:,2) , pf_idx{k,1}(:,3)); 
        %plot(pf_idx{k,1}(:,1), pf_idx{k,1}(:,2) , 'y-');
    end
    xlabel('c1');
    ylabel('c2'); 

   
    handles.pf_idx = pf_idx;
    handles.X = X;
       
    
    guidata(hObject, handles);
    
function tictoc2_Callback(hObject, eventdata, handles)
% hObject    handle to tictoc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tictoc2 as text
%        str2double(get(hObject,'String')) returns contents of tictoc2 as a double


% --- Executes during object creation, after setting all properties.
function tictoc2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tictoc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function QueryName3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to QueryName3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes during object creation, after setting all properties.
function axes12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes12


% --- Executes on selection change in hashCodeSelection_3d.
function hashCodeSelection_3d_Callback(hObject, eventdata, handles)
% hObject    handle to hashCodeSelection_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns hashCodeSelection_3d contents as cell array
%        contents{get(hObject,'Value')} returns selected item from hashCodeSelection_3d



video_dir = handles.video_dir;
data_dir  = handles.data_dir;

% load([data_dir '/filenames']); % File names
% load([data_dir '/targets']);   % Labels

hashCode_index = get(handles.hashCodeSelection_3d, 'Value');

switch hashCode_index
           
    case 1
        load([data_dir '/hashCodes_128']); 
        data = hashCodes_128;
        
    case 2
        load([data_dir '/hashCodes_256']); 
        data = hashCodes_256;
        
    case 3
        load([data_dir '/hashCodes_512']); 
        data = hashCodes_512;
       
    case 4
        load([data_dir '/hashCodes_1024']); 
        data = hashCodes_1024;
            
end


handles.data = data;
guidata(hObject, handles);






% --- Executes during object creation, after setting all properties.
function hashCodeSelection_3d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hashCodeSelection_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_pp_3d_Callback(hObject, eventdata, handles)
% hObject    handle to num_pp_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_pp_3d as text
%        str2double(get(hObject,'String')) returns contents of num_pp_3d as a double


% --- Executes during object creation, after setting all properties.
function num_pp_3d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_pp_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
