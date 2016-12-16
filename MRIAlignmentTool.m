function varargout = MRIAlignmentTool(varargin)
% MRIALIGNMENTTOOL MATLAB code for MRIAlignmentTool.fig
%      MRIALIGNMENTTOOL, by itself, creates a new MRIALIGNMENTTOOL or raises the existing
%      singleton*.
%
%      H = MRIALIGNMENTTOOL returns the handle to a new MRIALIGNMENTTOOL or the handle to
%      the existing singleton*.
%
%      MRIALIGNMENTTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MRIALIGNMENTTOOL.M with the given input arguments.
%
%      MRIALIGNMENTTOOL('Property','Value',...) creates a new MRIALIGNMENTTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MRIAlignmentTool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MRIAlignmentTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MRIAlignmentTool

% Last Modified by GUIDE v2.5 09-Dec-2016 16:35:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MRIAlignmentTool_OpeningFcn, ...
                   'gui_OutputFcn',  @MRIAlignmentTool_OutputFcn, ...
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


% --- Executes just before MRIAlignmentTool is made visible.
function MRIAlignmentTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MRIAlignmentTool (see VARARGIN)

% Choose default command line output for MRIAlignmentTool
handles.output = hObject;

setenv('FSLDIR','/usr/local/fsl');
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
% Update handles structure
axes(handles.axes1);
plot(0)
setAlignSection(hObject, 'off');
setModelSection(hObject, 'off');
setExportSection(hObject, 'off');
guidata(hObject, handles);

% UIWAIT makes MRIAlignmentTool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MRIAlignmentTool_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function [fullPath, slices] = getslicesFromPath(folderPath)
  currentDir = pwd;
  cd(folderPath)

  filecell = struct2cell([dir('*DCM*') dir('*MRDC*')]);

  fullPath = pwd;
  slices = filecell;
  cd(currentDir)

function checkEnableRearrange(handles)
  % check if struct has a specific field using isfield
  if isfield(handles,'upperPath') || isfield(handles,'lowerPath')
    set(handles.rearrangeButton, 'Enable', 'on')
  end

% --- Executes on button press in load_upper.
function load_upper_Callback(hObject, eventdata, handles)
  folder_path = uigetdir(pwd, 'Select a folder containing DICOM slices');
  if folder_path ~= 0
    handles = guidata(hObject);
    [upperPath, upperSlices] = getslicesFromPath(folder_path);
    handles.upperPath = upperPath;
    handles.upperSlices = upperSlices;

    set(handles.upperArmStatus, 'String', [num2str(length(handles.upperSlices)) ' DICOM slices detected'])
    checkEnableRearrange(handles);
    set(handles.load_lower, 'Enable', 'on');
    guidata(hObject, handles);
  end

% --- Executes on button press in load_lower.
function load_lower_Callback(hObject, eventdata, handles)
  folder_path = uigetdir(pwd, 'Select a folder containing DICOM slices');
  if folder_path ~= 0
    handles = guidata(hObject);
    [lowerPath, lowerSlices] = getslicesFromPath(folder_path);
    handles.lowerPath = lowerPath;
    handles.lowerSlices = lowerSlices;

    set(handles.lowerArmStatus, 'String', [num2str(length(handles.lowerSlices)) ' DICOM slices detected'])
    checkEnableRearrange(handles);
    guidata(hObject, handles);
  end

function indices = sortFileByName(files)
  if isempty(files)
    indices = [];
  else
    extCell = cellfun(@(x) x(end - 1),cellfun(@(x) strsplit(x,'.'),files(1,:),'UniformOutput',false));
    [~, indices] = sort(cellfun(@(x) str2num(x), extCell),'ascend');
  end

function indices = sortFileByExtension(files)
  if isempty(files)
    indices = [];
  else
    extCell = cellfun(@(x) x(end),cellfun(@(x) strsplit(x,'.'),files(1,:),'UniformOutput',false));
    [~, indices] = sort(cellfun(@(x) str2num(x), extCell),'ascend');
  end

% --- Executes on button press in rearrangeButton.
function rearrangeButton_Callback(hObject, eventdata, handles)
  % hObject    handle to rearrangeButton (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  handles = guidata(hObject);


  if ~isfield(handles, 'lowerSlices')
    choice = questdlg('Only one segment is loaded. Continue?', 'Warning', 'Cancel', 'Yes', 'Cancel');
    switch  choice
    case 'Cancel'
      return
    case 'Yes'
      handles.oneSegment = true;
      handles.lowerSlices = {};
    end
  else
    if isfield(handles, 'oneSegment')
      handles = rmfield(handles, 'oneSegment')
    end
  end

  if ~exist('wholeArm')
    mkdir wholeArm
  else
    cd wholeArm
    delete('*.MRDC');
    cd ..
  end

  setAlignSection(hObject, 'off');
  setModelSection(hObject, 'off');
  setExportSection(hObject, 'off');

  h = waitbar(0, 'Rearranging MRI slices...');% 'WindowStyle', 'modal');
  % set(handles.figure1, 'pointer', 'watch')

  sortExtension = get(handles.sortExtension, 'Value');

  sortedUpperIndex = zeros(1,1);
  sortedForeIndex = zeros(1,1);

  if sortExtension
    sortedUpperIndex = sortFileByExtension(handles.upperSlices);
    sortedForeIndex = sortFileByExtension(handles.lowerSlices);
  else
    sortedUpperIndex = sortFileByName(handles.upperSlices);
    sortedForeIndex = sortFileByName(handles.lowerSlices);
  end


  totalLength = length(handles.upperSlices) + length(handles.lowerSlices);
  for i = 1:length(handles.upperSlices)
    copyfile([handles.upperPath '/' char(handles.upperSlices(1, sortedUpperIndex(i)))], ['wholeArm/' sprintf('%04d.MRDC',i)])
    waitbar(i/totalLength);
  end

  for i = 1:length(handles.lowerSlices)
    copyfile([handles.lowerPath '/' char(handles.lowerSlices(1, sortedForeIndex(i)))], ['wholeArm/' sprintf('%04d.MRDC', i + length(handles.upperSlices))]);
    waitbar(i + length(handles.upperSlices)/totalLength);
  end

  [folderPath, allSlices] = getslicesFromPath('wholeArm');
  sortedIndex = sortFileByName(allSlices);

  fileNames = cell(length(allSlices),1);

  for i = 1:length(sortedIndex)
    fileNames{i} = char(strcat(folderPath, '/', allSlices(1, sortedIndex(i))));
  end

  handles.filePaths = fileNames;
  handles.alignmentOccured = false;

  axes(handles.axes1);
  % set(handles.figure1, 'CurrentAxes', handles.axes1);
  montage(fileNames,'DisplayRange',[])
  close(h);

  setAlignSection(hObject, 'on');

  set(handles.upperframeSlider, 'Min', 1);
  set(handles.upperframeSlider, 'Max', length(handles.upperSlices));
  set(handles.upperframeSlider, 'Value', length(handles.upperSlices));
  set(handles.upperNumberField, 'String', length(handles.upperSlices));
  % isfield(handles, 'oneSegment')
  if ~isfield(handles, 'oneSegment')
    set(handles.lowerframeSlider, 'Min', length(handles.upperSlices) + 1);
    set(handles.lowerframeSlider, 'Max', totalLength);
    set(handles.lowerframeSlider, 'Value', length(handles.upperSlices) + 1);
    set(handles.lowerNumberField, 'String', num2str(length(handles.upperSlices) + 1));
  else
    set(handles.lowerframeSlider, 'Enable', 'off');
    set(handles.lowerNumberField, 'Enable', 'off');
  end

  set(handles.cpalignButton, 'Enable', 'off');
  set(handles.autoalignButton, 'Enable', 'off');
  % set(handles.mergeButton, 'Enable', 'off');

  guidata(hObject, handles);

  % set(handles.figure1, 'pointer', 'arrow')

function upperNumberField_Callback(hObject, eventdata, handles)
  display2slices(hObject);
  currentVal = get(hObject, 'String');
  set(handles.upperframeSlider, 'Value', str2num(currentVal));
function lowerNumberField_Callback(hObject, eventdata, handles)
  display2slices(hObject);
  currentVal = get(hObject, 'String');
  set(handles.lowerframeSlider, 'Value', str2num(currentVal));

% --- Executes on slider movement.
function upperframeSlider_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'Value');
  set(handles.upperNumberField, 'String', num2str(round(currentVal)));
  display2slices(hObject);
% --- Executes on slider movement.
function lowerframeSlider_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'Value');
  set(handles.lowerNumberField, 'String', num2str(round(currentVal)));
  display2slices(hObject);

% --- Executes on button press in view2sliceButton.
function view2sliceButton_Callback(hObject, eventdata, handles)
  display2slices(hObject);


function display2slices(hObject)
  % --- display 2 slices together
  handles = guidata(hObject);
  upperNum = str2num(get(handles.upperNumberField, 'String'));
  lowerNum = str2num(get(handles.lowerNumberField, 'String'));

  assignin('base','filePaths',handles.filePaths);
  assignin('base','upperNum',upperNum);
  assignin('base','lowerNum',lowerNum);

  base_upperInfo = dicominfo(handles.filePaths{upperNum});
  [base_upperImg, ~] = dicomread(base_upperInfo);
  base_upperImg = mat2gray(base_upperImg);
  handles.upperImg = base_upperImg;

  set(handles.axes1, 'Visible', 'on');
  axes(handles.axes1)

  % upperImg = imadjust(base_upperImg,[0 20/255],[0 1]);
  if ~isfield(handles, 'oneSegment')
    base_foreInfo = dicominfo(handles.filePaths{lowerNum});
    [base_foreImg, ~] = dicomread(base_foreInfo);
    base_foreImg = mat2gray(base_foreImg);
    % foreImg = imadjust(base_foreImg,[0 20/255],[0 1]);
    handles.foreImg = base_foreImg;
    imshowpair(base_upperImg, base_foreImg, 'montage');
    set(handles.cpalignButton, 'Enable', 'on');
    % set(handles.mergeButton, 'Enable', 'on');
  else
    imshow(base_upperImg);
  end

  guidata(hObject, handles);


% --- Executes on button press in cpalignButton.
function cpalignButton_Callback(hObject, eventdata, handles)
  handles = guidata(hObject);

  [movingout, fixedout] = cpselect(handles.foreImg, handles.upperImg, 'Wait', true);
  tform = fitgeotrans(movingout,fixedout, 'affine');

  Roriginal = imref2d(size(handles.upperImg));
  recovered = imwarp(handles.foreImg, tform,'OutputView',Roriginal);

  axes(handles.axes1);
  imshowpair(handles.upperImg, recovered);

  handles.transform2d = tform;
  handles.alignmentOccured = true;
  % set(handles.mergeButton, 'Enable', 'on');

  guidata(hObject, handles);

% --- Executes on button press in autoalignButton.
function autoalignButton_Callback(hObject, eventdata, handles)
% hObject    handle to autoalignButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in mergeButton.
function mergeButton_Callback(hObject, eventdata, handles)
% hObject    handle to mergeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  handles = guidata(hObject);
  h = waitbar(0, 'Merging MRI slices...');

  totalLength = length(handles.filePaths);
  upperIndex = str2num(get(handles.upperNumberField, 'String'));
  foreIndex = str2num(get(handles.lowerNumberField, 'String'));

  if isfield(handles, 'oneSegment')
    upperIndex = totalLength;
  else
    display2slices(hObject);
    handles = guidata(hObject);
    Roriginal = imref2d(size(handles.upperImg));
  end

  for i = 1:upperIndex
    upperInfo = dicominfo(handles.filePaths{i});
    upperImg = dicomread(upperInfo);
    upperInfo.SeriesDescription = ['vol_slice_' sprintf('%0.4d', i)];
    dicomwrite(upperImg,'temp_vol', upperInfo);
    dicm2nii('temp_vol',['wholeArmNii'],'nii');
    waitbar(i/totalLength);
  end
  delete('temp_vol');

  if ~isfield(handles, 'oneSegment')
    foreFileIndex = upperIndex + 1;
    % transform the rest of the forearm images, convert to nii and move to another folder
    for i = foreIndex:length(handles.filePaths)
      foreInfo = dicominfo(handles.filePaths{i});
      [foreImg, foreMap] = dicomread(foreInfo);

      if handles.alignmentOccured == true,
        recovered = imwarp(foreImg, handles.transform2d,'OutputView',Roriginal);
      else
        recovered = foreImg;
      end
      foreInfo.SeriesDescription = ['vol_slice_' sprintf('%0.4d', foreFileIndex)];
      dicomwrite(recovered, 'transformed_vol', foreInfo);
      dicm2nii('transformed_vol',['wholeArmNii'],'nii');
      foreFileIndex = foreFileIndex + 1;

      waitbar(i/totalLength);
    end
    % delete temp filecell
    delete('transformed_vol');
  end

  curdatetime = fix(clock);
  curdatetime = arrayfun(@(x) num2str(x), curdatetime, 'UniformOutput', false);
  curdatetime = strjoin(curdatetime(1:5), '_');
  handles.mergedFileName = ['wholeArm_' curdatetime];

  cd('wholeArmNii')
  system(char(['/usr/local/fsl/bin/fslmerge -z ' handles.mergedFileName ' vol_slice_*']))
  gunzip([handles.mergedFileName '.nii.gz'])
  numSlices = length(dir('vol_slice_*'));
  delete('vol_slice_*');
  set(handles.load3DModel,'enable','on');
  cd ..

  close(h);
  guidata(hObject, handles);
  load3DModel_Callback(hObject, [], handles);

% --- Executes on button press in load3DModel.
function load3DModel_Callback(hObject, eventdata, handles)
% hObject    handle to load3DModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

set(handles.figure1, 'pointer', 'watch')

cd('wholeArmNii')
face_param = struct;
face_param.radius = [-2 2 3 -3 -3];
face_parm.step = 1;
face_parm.pmax = 0.099;

[Fface, Vface] = face_extract([handles.mergedFileName '.nii']);
% imagefile = 'wholeArmNii/Ax_T2_Armtest.nii';
% [Fface, Vface] = face_extract('wholeArmNii/wholeArm_2016_12_9_15_54.nii');


surf_face.V = Vface * 1000; %vbmeg converts mm to m, convert back
surf_face.F = Fface;
surf_face.face_parm = face_parm;

Nmri = 5000;
[surf_face.F_reduce, surf_face.V_reduce] = reducepatch(surf_face.F, surf_face.V, 2 * Nmri);
axes(handles.axes1);

handles.modelPlot = trisurf(surf_face.F_reduce, surf_face.V_reduce(:,1), surf_face.V_reduce(:,2), surf_face.V_reduce(:,3), 'FaceAlpha', 0.1, 'EdgeColor', [0.5 0.5 0.5])
% handles.modelPlot = trisurf(surf_face.F, surf_face.V(:,1), surf_face.V(:,2), surf_face.V(:,3), 'FaceAlpha', 0.1, 'EdgeColor', [0.5 0.5 0.5])

axis equal
set( gcf, 'menubar', 'figure' )
rotate3d on

handles.surf_face = surf_face;
guidata(hObject, handles);
set(handles.figure1, 'pointer', 'arrow')
set(handles.loadsensorButton, 'Enable', 'on');
cd ..

% --- Executes on button press in loadsensorButton.
function loadsensorButton_Callback(hObject, eventdata, handles)
  % hObject    handle to loadsensorButton (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  handles = guidata(hObject);
  [file_name, file_path] = uigetfile('*.csv', 'Select a sensor file');

  if file_name ~= 0
    [labels, pos] = readOptoPositions([file_path file_name]);
    electrodes = pos(:,:);

    raw_electrodes = electrodes;

    temp_electrodes = electrodes;
    electrodes(2,:) = temp_electrodes(7,:);
    electrodes(3,:) = temp_electrodes(2,:);
    electrodes(4,:) = temp_electrodes(3,:);
    electrodes(6,:) = temp_electrodes(4,:);
    electrodes(7,:) = temp_electrodes(6,:);

    handles.raw_electrodes = raw_electrodes;
    handles.electrodes = electrodes;
    handles.labels = labels;

    hold on
    handles.hRef = scatter3(electrodes(:,1), electrodes(:,2), electrodes(:,3),'filled','lineWidth',20);
    handles.hE = scatter3(electrodes(1:7,1),electrodes(1:7,2),electrodes(1:7,3), 'filled','lineWidth',20);
    handles.hT = text(electrodes(:,1),electrodes(:,2),electrodes(:,3), num2str(handles.labels));
    % update handle object
    hold off

    % msgbox('Sensor positions loaded');

    set(handles.markmodeButton, 'Enable', 'on');
  end

  handles.rotateMat = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
  handles.scaleMat = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
  handles.translateMat = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];

  guidata(hObject, handles);

  % --- Executes on button press in markmodeButton.
  function markmodeButton_Callback(hObject, eventdata, handles)
  % hObject    handle to markmodeButton (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  % Hint: get(hObject,'Value') returns toggle state of markmodeButton
  handles = guidata(hObject);
  set( gcf, 'menubar', 'figure' )
  dcm_obj = datacursormode(gcf);
  set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on');
  handles.dcm_obj = dcm_obj;
  set(handles.transformButton, 'Enable', 'on');
guidata(hObject, handles);

% --- Executes on button press in transformButton.
function transformButton_Callback(hObject, eventdata, handles)
% hObject    handle to transformButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
set(handles.figure1, 'pointer', 'watch')

c_info = getCursorInfo(handles.dcm_obj)
assignin('base','c_info',c_info);

a_markers = zeros(length(c_info), 3);
index = fliplr(1:length(c_info));
for i = 1:length(c_info)
  a_markers(i,1:3) = c_info(index(i)).Position;
end
assignin('base','index', index);
assignin('base','a_markers',a_markers);
assignin('base','electrodes',handles.electrodes);

te_electrodes = cpd_align(handles.electrodes, a_markers);
% update plot
axes(handles.axes1);

delete(handles.hRef);
delete(handles.hE);
delete(handles.hT);

hold on
handles.hRef = scatter3(te_electrodes(:,1), te_electrodes(:,2), te_electrodes(:,3),'filled','lineWidth',20);
handles.hE = scatter3(te_electrodes(1:7,1),te_electrodes(1:7,2),te_electrodes(1:7,3), 'filled','lineWidth',20);
handles.hT = text(te_electrodes(:,1),te_electrodes(:,2),te_electrodes(:,3), num2str(handles.labels));
% update handle object
hold off
handles.electrodes = te_electrodes;

setModelSection(hObject, 'on');
guidata(hObject, handles);
set(handles.figure1, 'pointer', 'arrow')

function transformed_electrodes = cpd_align(electrodes, a_markers)
% Set the options
opt.method='rigid'; % use rigid registration
opt.viz=0;          % show every iteration
opt.outliers=0;     % do not assume any noise

opt.normalize=0;    % normalize to unit variance and zero mean before registering (default)
opt.scale=1;        % estimate global scaling too (default)
opt.rot=1;          % estimate strictly rotational matrix (default)
opt.corresp=0;      % do not compute the correspondence vector at the end of registration (default)

opt.max_it=100;     % max number of iterations
opt.tol=1e-8;       % tolerance

[rows, ~] = size(a_markers);
e_markers = electrodes(1:rows, 1:3);

[ret_r,ret_t] = rigid_transform_3D(e_markers, a_markers);

te_markers = (ret_r*e_markers') + repmat(ret_t, 1, rows);
te_markers = te_markers';
te_electrodes = (ret_r * electrodes') + repmat(ret_t, 1, length(electrodes));
te_electrodes = te_electrodes';

[Transform, C]=cpd_register(a_markers, te_electrodes, opt);
% transformed_electrodes = Transform.Y;
transformed_electrodes = te_electrodes;

function rotate(hObject, about, degree)
  handles = guidata(hObject);
  tMat = zeros(4,4);
  rotRad = degree * pi / 180;

  if isfield(handles, 'prevRMat') == false
    handles.prevRMat = [0 0 0];
  end

  switch about
    case 'x'
      difRot = rotRad - handles.prevRMat(1);
      tMat = [1 0 0 0; 0 cos(difRot) sin(difRot) 0; 0 -sin(difRot) cos(difRot) 0; 0 0 0 1];
      handles.prevRMat(1) = rotRad;
    case 'y'
      difRot = rotRad - handles.prevRMat(2);
      tMat = [cos(difRot) 0 -sin(difRot) 0; 0 1 0 0; sin(difRot) 0 cos(difRot) 0; 0 0 0 1];
      handles.prevRMat(2) = rotRad;
    case 'z'
      difRot = rotRad - handles.prevRMat(3);
      tMat = [cos(difRot) sin(difRot) 0 0; -sin(difRot) cos(difRot) 0 0; 0 0 1 0; 0 0 0 1];
      handles.prevRMat(3) = rotRad;
    otherwise
      tMat = ones(4,4);
  end

  handles.rotateMat = tMat;
  electrodes = handles.electrodes;
  electrodes = [electrodes ones(length(electrodes), 1)] * tMat;
  handles.electrodes = electrodes(:,1:3);
  te_electrodes = electrodes * handles.translateMat * handles.scaleMat;

  delete(handles.hRef);
  delete(handles.hE);
  delete(handles.hT);
  hold on
  handles.hRef = scatter3(te_electrodes(:,1), te_electrodes(:,2), te_electrodes(:,3),'filled','lineWidth',20);
  handles.hE = scatter3(te_electrodes(1:7,1),te_electrodes(1:7,2),te_electrodes(1:7,3), 'filled','lineWidth',20);
  handles.hT = text(te_electrodes(:,1),te_electrodes(:,2),te_electrodes(:,3), num2str(handles.labels));
  hold off
  % comment out for beautiful effects
  guidata(hObject, handles);

% --- Executes on slider movement.
function rotateX_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'Value');
  set(handles.rotateXField, 'String', num2str(currentVal));
  rotate(hObject, 'x', currentVal);
% --- Executes on slider movement.
function rotateY_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'Value');
  set(handles.rotateYField, 'String', num2str(currentVal));
  rotate(hObject, 'y', currentVal);

% --- Executes on slider movement.
function rotateZ_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'Value');
  set(handles.rotateZField, 'String', num2str(currentVal));
  rotate(hObject, 'z', currentVal);

function rotateXField_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'String');
  rotate(hObject, 'x', str2num(currentVal));
  set(handles.rotateX, 'Value', str2num(currentVal));

function rotateYField_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'String');
  rotate(hObject, 'y', str2num(currentVal));
  set(handles.rotateY, 'Value', str2num(currentVal));

function rotateZField_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'String');
  rotate(hObject, 'z', str2num(currentVal));
  set(handles.rotateZ, 'Value', str2num(currentVal));

% TODO: incremental translation + rotation

function translate(hObject, along, matrix)
  handles = guidata(hObject);
  electrodes = handles.electrodes;

  tMat = handles.translateMat;
  switch along
    case 'x'
      tMat(4,1) = matrix(1);
    case 'y'
      tMat(4,2) = matrix(2);
    case 'z'
      tMat(4,3) = matrix(3);
  end

  % tMat = [1 0 0 0; 0 1 0 0; 0 0 1 0; difMat(1) difMat(2) difMat(3) 1];
  te_electrodes = [electrodes ones(length(electrodes), 1)] * tMat * handles.scaleMat;
  handles.translateMat = tMat;

  delete(handles.hRef);
  delete(handles.hE);
  delete(handles.hT);
  hold on
  handles.hRef = scatter3(te_electrodes(:,1), te_electrodes(:,2), te_electrodes(:,3),'filled','lineWidth',20);
  handles.hE = scatter3(te_electrodes(1:7,1),te_electrodes(1:7,2),te_electrodes(1:7,3), 'filled','lineWidth',20);
  handles.hT = text(te_electrodes(:,1),te_electrodes(:,2),te_electrodes(:,3), num2str(handles.labels));
  hold off

  % comment out for beautiful effects
  guidata(hObject, handles);

% --- Executes on slider movement.
function translateX_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'Value');
  set(handles.translateXField, 'String', num2str(currentVal));
  translate(hObject, 'x', [currentVal 0 0]);

% --- Executes on slider movement.
function translateY_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'Value');
  set(handles.translateYField, 'String', num2str(currentVal));
  translate(hObject, 'y', [0 currentVal 0]);

% --- Executes on slider movement.
function translateZ_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'Value');
  set(handles.translateZField, 'String', num2str(currentVal));
  translate(hObject, 'z', [0 0 currentVal]);

function translateXField_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'String');
  translate(hObject, 'x', [str2num(currentVal) 0 0]);
  set(handles.translateX, 'Value', str2num(currentVal));

function translateYField_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'String');
  translate(hObject, 'y', [0 str2num(currentVal) 0]);
  set(handles.translateX, 'Value', str2num(currentVal));

function translateZField_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'String');
  translate(hObject, 'z', [0 0 str2num(currentVal)]);
  set(handles.translateX, 'Value', str2num(currentVal));

function scaler(hObject, along, matrix)
  handles = guidata(hObject);
  electrodes = handles.electrodes;

  tMat = handles.scaleMat;
  switch along
    case 'x'
      tMat(1,1) = matrix(1);
    case 'y'
      tMat(2,2) = matrix(2);
    case 'z'
      tMat(3,3) = matrix(3);
    otherwise
      tMat(1,1) = matrix(1);
      tMat(2,2) = matrix(2);
      tMat(3,3) = matrix(3);
  end

  % tMat = [difMat(1) 0 0 0; 0 difMat(2) 0 0; 0 0 difMat(3) 0; 0 0 0 1];
  te_electrodes = [electrodes ones(length(electrodes), 1)] * handles.translateMat * tMat;
  handles.scaleMat = tMat;

  delete(handles.hRef);
  delete(handles.hE);
  delete(handles.hT);
  hold on
  handles.hRef = scatter3(te_electrodes(:,1), te_electrodes(:,2), te_electrodes(:,3),'filled','lineWidth',20);
  handles.hE = scatter3(te_electrodes(1:7,1),te_electrodes(1:7,2),te_electrodes(1:7,3), 'filled','lineWidth',20);
  handles.hT = text(te_electrodes(:,1),te_electrodes(:,2),te_electrodes(:,3), num2str(handles.labels));
  hold off

  % comment out for beautiful effects
  guidata(hObject, handles);


% TODO: update x, y, & z accordingly when updating the scale as a whole
% --- Executes on slider movement.
function scaleAll_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'Value');
  set(handles.scaleAllField, 'String', num2str(currentVal));

  set(handles.scaleXField, 'String', num2str(currentVal));
  set(handles.scaleYField, 'String', num2str(currentVal));
  set(handles.scaleZField, 'String', num2str(currentVal));

  set(handles.scaleX, 'Value', currentVal);
  set(handles.scaleY, 'Value', currentVal);
  set(handles.scaleZ, 'Value', currentVal);

  scaler(hObject, 'all', [currentVal currentVal currentVal]);

% --- Executes on slider movement.
function scaleX_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'Value');
  set(handles.scaleXField, 'String', num2str(currentVal));
  scaler(hObject, 'x', [currentVal 0 0]);

% --- Executes on slider movement.
function scaleY_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'Value');
  set(handles.scaleYField, 'String', num2str(currentVal));
  scaler(hObject, 'y', [0 currentVal 0]);

  % --- Executes on slider movement.
  function scaleZ_Callback(hObject, eventdata, handles)
    currentVal = get(hObject, 'Value');
    set(handles.scaleZField, 'String', num2str(currentVal));
    scaler(hObject, 'z', [0 0 currentVal]);


function scaleAllField_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'String');
  currentVal = str2num(currentVal);
  scaler(hObject, 'all', [currentVal currentVal currentVal]);
  set(handles.scaleAll, 'Value', currentVal);

function scaleXField_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'String');
  currentVal = str2num(currentVal);
  scaler(hObject, 'x', [currentVal 0 0]);
  set(handles.scaleX, 'Value', currentVal);

function scaleYField_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'String');
  currentVal = str2num(currentVal);
  scaler(hObject, 'y', [0 currentVal 0]);
  set(handles.scaleY, 'Value', currentVal);

function scaleZField_Callback(hObject, eventdata, handles)
  currentVal = get(hObject, 'String');
  currentVal = str2num(currentVal);
  scaler(hObject, 'z', [0 0 currentVal]);
  set(handles.scaleZ, 'Value', currentVal);


% --- Executes on button press in wrapSensorPosition.
function wrapSensorPosition_Callback(hObject, eventdata, handles)
% hObject    handle to wrapSensorPosition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  handles = guidata(hObject);
  electrodes = handles.electrodes;

  h = waitbar(0, 'Positioning sensor locations');

  te_electrodes = [electrodes ones(length(electrodes), 1)] * handles.translateMat * handles.scaleMat;
  te_electrodes = te_electrodes(:,1:3);

  wrappedSensor = wrapSensor(handles.surf_face.V, te_electrodes);

  close(h);

  delete(handles.hRef);
  delete(handles.hE);
  delete(handles.hT);
  hold on

  expand = 1;
  % tMat = [expand 0 0 0; 0 expand 0 0; 0 0 expand 0; 0 0 0 1];
  wrappedSensor = [wrappedSensor ones(length(wrappedSensor), 1)];% * tMat;

  handles.hRef = scatter3(wrappedSensor(:,1), wrappedSensor(:,2), wrappedSensor(:,3),'filled','lineWidth',20);
  handles.hE = scatter3(wrappedSensor(1:7,1),wrappedSensor(1:7,2),wrappedSensor(1:7,3), 'filled','lineWidth',20);

  handles.wrapped_electrodes = wrappedSensor;

  wrappedSensor = wrappedSensor(:,1:3);
  expand = 1.05;
  tMat = [expand 0 0 0; 0 expand 0 0; 0 0 expand 0; 0 0 0 1];
  wrappedSensor = [wrappedSensor ones(length(wrappedSensor), 1)] * tMat;

  handles.hT = text(wrappedSensor(:,1),wrappedSensor(:,2),wrappedSensor(:,3), num2str(handles.labels));
  hold off

  setExportSection(hObject, 'on');
  guidata(hObject, handles);

function wrappedSensor = wrapSensor(meshData, electrodeData)
  [rows, ~] = size(electrodeData);

  wrappedSensor = zeros(rows, 3);
  for i=1:rows
    distances = sqrt(sum( bsxfun(@minus, meshData, electrodeData(i,:)).^2, 2));
    [r, c, val] = find(distances==min(distances));
    wrappedSensor(i,:) = meshData(r(1),:);
  end

function setAlignSection(hObject, status)
  handles = guidata(hObject);
  set(handles.upperNumberField, 'Enable', status);
  set(handles.upperframeSlider, 'Enable', status);
  set(handles.lowerNumberField, 'Enable', status);
  set(handles.lowerframeSlider, 'Enable', status);
  set(handles.view2sliceButton, 'Enable', status);
  set(handles.cpalignButton, 'Enable', status);
  set(handles.mergeButton, 'Enable', status);
  guidata(hObject, handles);

function setModelSection(hObject, status)
  handles = guidata(hObject);
  set(handles.load3DModel, 'Enable', status);
  set(handles.loadsensorButton, 'Enable', status);
  set(handles.markmodeButton, 'Enable', status);
  set(handles.transformButton, 'Enable', status);
  set(handles.rotateX, 'Enable', status);
  set(handles.rotateY, 'Enable', status);
  set(handles.rotateZ, 'Enable', status);
  set(handles.rotateXField, 'Enable', status);
  set(handles.rotateYField, 'Enable', status);
  set(handles.rotateZField, 'Enable', status);
  set(handles.translateX, 'Enable', status);
  set(handles.translateY, 'Enable', status);
  set(handles.translateZ, 'Enable', status);
  set(handles.translateXField, 'Enable', status);
  set(handles.translateYField, 'Enable', status);
  set(handles.translateZField, 'Enable', status);
  set(handles.scaleAll, 'Enable', status);
  set(handles.scaleX, 'Enable', status);
  set(handles.scaleY, 'Enable', status);
  set(handles.scaleZ, 'Enable', status);
  set(handles.scaleAllField, 'Enable', status);
  set(handles.scaleXField, 'Enable', status);
  set(handles.scaleYField, 'Enable', status);
  set(handles.scaleZField, 'Enable', status);
  set(handles.wrapSensorPosition, 'Enable', status);
  guidata(hObject, handles);

function setExportSection(hObject, status)
  handles = guidata(hObject);
  set(handles.export_data, 'Enable', status);
  set(handles.saveSensorPosition, 'Enable', status);
  set(handles.savemodelButton, 'Enable', status);
  set(handles.saveFigure, 'Enable', status);
  guidata(hObject, handles);

  % --- Executes on button press in savemodelButton.
  function savemodelButton_Callback(hObject, eventdata, handles)
  % hObject    handle to savemodelButton (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
   handles = guidata(hObject);
   [fileName, pathName] = uiputfile(['arm_model_surface_' strrep(date,'-','_') '.mat'],'Save file name');
   currentDirectory = pwd;
   cd(pathName);
   fileName
   if isfield(handles, 'importPlot')
     model = struct;
     model.skin = handles.surf_face;
     model.bone = handles.bone_face;
     model.muscle = handles.muscle_face;
   else
     model = handles.model;
   end

   save(fileName, 'model');
   cd(currentDirectory);

  % --- Executes on button press in saveSensorPosition.
  function saveSensorPosition_Callback(hObject, eventdata, handles)
  % hObject    handle to saveSensorPosition (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  handles = guidata(hObject);
  [fileName, pathName] = uiputfile(['sensor_position_' strrep(date,'-','_') '.csv'],'Save file name');
  currentDirectory = pwd;
  cd(pathName);
  csvwrite(fileName,[handles.labels handles.wrapped_electrodes(:,1:3)]);
  cd(currentDirectory);


  % --- Executes on button press in saveFigure.
  function saveFigure_Callback(hObject, eventdata, handles)
  % hObject    handle to saveFigure (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  % https://www.mathworks.com/matlabcentral/answers/86693-save-axes-plot-as-fig-in-a-gui

  Fig2 = figure('visible','off');
  copyobj(handles.axes1, Fig2);
  savefig(Fig2, 'myFigure.fig')

  handles = guidata(hObject);
  [fileName, pathName] = uiputfile(['arm_model_' strrep(date,'-','_') '.figure'],'Save file name');
  currentDirectory = pwd;
  cd(pathName);
  fileName
  savefig(gcf, fileName);
  cd(currentDirectory);

  % --- Executes on button press in export_data.
  function export_data_Callback(hObject, eventdata, handles)
  % hObject    handle to export_data (see GCBO)
  % eventdata  reserved - to be defined in a future version of MATLAB
  % handles    structure with handles and user data (see GUIDATA)
  handles = guidata(hObject);
  [fileName, pathName] = uiputfile(['arm_data_' strrep(date,'-','_') '.mat'],'Save file name')
  currentDirectory = pwd;
  cd(pathName);
  export_data = struct;
  export_data.surf_face = handles.surf_face;
  export_data.electrodes = handles.wrapped_electrodes(:,1:3);
  export_data.raw_electrodes = handles.raw_electrodes(:,1:3);
  export_data.labels = handles.labels;
  save(fileName, 'export_data');
  cd(currentDirectory);

% --------------------------------------------------------------------
function About_menu_Callback(hObject, eventdata, handles)
% hObject    handle to About_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Help_menu_Callback(hObject, eventdata, handles)
% hObject    handle to Help_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in makeModelTransparent.
function makeModelTransparent_Callback(hObject, eventdata, handles)
% hObject    handle to makeModelTransparent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
hold on
surf_face = handles.surf_face;

if isfield(handles, 'modelPlot')
  delete(handles.modelPlot);
  handles.modelPlot = trisurf(surf_face.F_reduce, surf_face.V_reduce(:,1), surf_face.V_reduce(:,2), surf_face.V_reduce(:,3), 'FaceAlpha', 0.1, 'EdgeColor', [0.5 0.5 0.5])
elseif isfield(handles, 'importPlot')
  delete(handles.importPlot);
  handles.importPlot = trisurf(surf_face.F, surf_face.V(:,1), surf_face.V(:,2), surf_face.V(:,3),  'FaceAlpha', 0.1, 'EdgeColor', [0.5 0.5 0.5])
end
axis equal
hold off
guidata(hObject, handles);

% --- Executes on button press in makeModelSkin.
function makeModelSkin_Callback(hObject, eventdata, handles)
% hObject    handle to makeModelSkin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
hold on
surf_face = handles.surf_face;

if isfield(handles, 'modelPlot')
  delete(handles.modelPlot);
  handles.modelPlot = trisurf(surf_face.F_reduce, surf_face.V_reduce(:,1), surf_face.V_reduce(:,2), surf_face.V_reduce(:,3), 'FaceAlpha', 1, 'FaceColor', [255 218 200]/255, 'EdgeColor', [255 218 200]/255, 'facelighting','gouraud', 'BackFaceLighting', 'unlit')
  delete(findall(gcf,'Type','light'))
  lightangle(-45,30)
  lightangle(-135,30)
  lightangle(45,-30)
  lightangle(135,-30)
elseif isfield(handles, 'importPlot')
  delete(handles.importPlot);
  handles.importPlot = trisurf(surf_face.F, surf_face.V(:,1), surf_face.V(:,2), surf_face.V(:,3), 'FaceAlpha', 0.9, 'FaceColor', [255 218 200]/255, 'EdgeColor', [248./255.,216./255.,183./255.], 'facelighting','gouraud', 'BackFaceLighting', 'lit')
end
axis equal
hold off


guidata(hObject, handles);


% --- Executes on button press in importMesh.
function importMesh_Callback(hObject, eventdata, handles)
% hObject    handle to importMesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file_name, file_path] = uigetfile('*.mat', 'Select a .mat file');

if file_name ~= 0,
  handles = guidata(hObject);
  load([file_path file_name]);

  axes(handles.axes1);
  % Plot together
  D = full_mask{1};
  data_size = size(D);

  if isfield(handles,'modelPlot')
    delete(modelPlot)
  end

  handles.surf_face.F = fskin(:,[1 2 3]);
  handles.surf_face.V = nskin(:,[1 2 3]);
  handles.bone_face.F = fbone(:,[1 2 3]);
  handles.bone_face.V = nbone(:,[1 2 3]);
  handles.muscle_face.F = fmuscle(:,[1 2 3]);
  handles.muscle_face.V = nmuscle(:,[1 2 3]);

  surf_face = handles.surf_face;

  hold on
  plotmesh(nbone(:,[1 2 3]),fbone(:,[1 2 3]),'facealpha',0.7,'facecolor', [1 1 1], 'EdgeAlpha', 0.25);
  plotmesh(nmuscle(:,[1 2 3]),fmuscle(:,[1 2 3]),'facealpha',0.5,'facecolor', [0.85 0 0], 'EdgeAlpha', 0.25);
  handles.importPlot = trisurf(surf_face.F, surf_face.V(:,1), surf_face.V(:,2), surf_face.V(:,3), 'FaceAlpha', 0.5, 'FaceColor', [255 218 200]/255, 'EdgeColor', [248./255.,216./255.,183./255.], 'facelighting','gouraud', 'BackFaceLighting', 'lit')

  % plotmesh(nskin(:,[1 2 3]),fskin(:,[1 2 3]),'facealpha',0.7,'facecolor', [248./255.,216./255.,183./255.], 'EdgeAlpha', 0.25);

  axis equal
  hold off


  set( gcf, 'menubar', 'figure' )
  guidata(hObject, handles);
end
% plotmesh(nskin(:,[1 2 3]),fskin(:,[1 2 3]),'facealpha',0.7,'facecolor', [248./255.,216./255.,183./255.], 'EdgeAlpha', 0.25);
% hold off
% legend('Bone','Muscle','Fat and Skin')
% xlabel('X [mm]')
% ylabel('Y [mm]')
% zlabel('Z [mm]')
% view(37.5,30)
% set(gca, 'fontsize',14)
% set(gca,'xtick',([0 40 60 80 100]))
% set(gcf,'color','w')
