classdef Luminosa_smMIET_parameters < handle
% Graphical user interface to input the sample parameters. This function is
% called by the main function "MIET_gui.m" when an image is to be evaluated 
% in single-molecule- or region-of-interest-mode.
%
% (c) Daja Ruhlandt 2014
    
    properties
        Figure          % "canvas" for placing all other objects
        PanelDispersion % "canvas" for deciding if calculation is monochromatic or uses spectrum
        PanelUpper      % "canvas" for placing everything concerning upper stack of materials
        PanelMolecule   % "canvas" for placing everything concerning the molecule's layer
        PanelLower      % "canvas" for placing everything concerning lower stack of materials
        PanelGeneral    % "canvas" for placing everything concerning general information
        
        DispText        % text describing difference between mono- & polychromatic
        DispChoose      % drop-down menu to choose: mono- or polychromatic?
        PolyInfotext    % text listing dummy refractive indices for polychromatic evaluation
        RefrIndexUpper  % editable table for entering indices of refraction (upper stack)
        RefrIndexLower  % editable table for entering indices of refraction (lower stack)
        ThicknessUpper  % editable table for entering thicknesses of layers (upper stack)
        ThicknessLower  % editable table for entering thicknesses of layers (lower stack)
        RefrIndexMolec  % editable text field for entering refractive index of molecule's layer
        ThicknessMolec  % editable text field for entering thickness of molecule's layer
        QuantumYield    % editable text field for free space QY
        DipoleOrientation % popup menu for the dipole-orientation model used in MIET calibration
        LaserMode       % drop-down menu to choose laser polarization mode
        NA              % editable text fied for entering numerical aperture of the objective
        Focus           % editable text field for entering defocussing [nm]
        ExcWavelength   % editable text field for entering excitation wavelength
        EmiWavelength   % editable text field for entering emission wavelength
        Wavel_Small     % smallest wavelength used in polychromatic evaluation (in nm)
        Wavel_Large     % largest wavelength used in polychromatic evaluation (in nm)
        LoadSpectrumBt  % button to open a text file containing the emission spectrum of the melocule
        SaveSettingsBt  % save refr. ind., thickness, z-parameters
        LoadSettingsBt  % load refr. ind., thickness, z-parameters
        ChooseMetal     % drop-down menu to choose metal for finding refractive index
        ReadRefrInd     % text field displaying the found refractive index
        ChooseCurveType % drop-down menu to choose: calibration curve to 1st maximum, 1st minumum or manual
        
        OKButton        % pushbutton to finish setting the parameters
        CancelButton    % pushbutton to cancel setting the parameters
        
        Label_n1        % labels for indices of refraction of upper stack
        Label_n1T
        Label_n1B
        Label_d1        % label for thicknesses of layers of upper stack
        Label_n         % label for index of refraction of molecule's layer
        Label_d         % label for thickness of molecule's layer
        Label_n0        % labels for indices of refraction of lower stack
        Label_n0T
        Label_n0B
        Label_d0        % label for thicknesses of layers of lower stack
        Label_laserMode % label for polarization of the laser (linear,radial,azimutal)
        Label_NA        % labels for numerical aperture of microscope objective
        Label_Focus     % label for defocussing 
        Label_ExcWavelen% label for excitation wavelength
        Label_EmiWavelen% label for emission wavelength
        Label_WavelDisp % labels for allowed wavelength range for polychromatic evaluation
        Label_WavelRange
        Label_SpectrumFile% name of the file containing the spectrum
        Label_FindRefrInd % label for finding refr. index of a material at certain wavelength
        Label_QY        % Label for free space qunatum yield
        Label_DipoleOrientation % label for the dipole-orientation model
        Label_ChooseCurveType
        
        Metals          % structure containing metal refractive indices
        heightIncr      % auxiliary variable for layout of app
        widthIncr       % auxiliary variable for layout of app
        SpectrumFile    % name of the file containing the spectrum
        Monochrome      % boolean: true=monochromatic evaluation, false=polychromatic evaluation
    end 

    
    methods
        % %  functions for designing and controlling the GUI  % % % % % % %
        function app = Luminosa_smMIET_parameters(varargin)  % "constructor" (can accept structure with parameter values as input)
        % auxiliary constants (change if window size is changed, see resizeApp()) 
            app.heightIncr = 34;
            app.widthIncr = 36;
            LayoutMolec = (3.6*app.heightIncr-120)/7;
        % 'invisible' variables
            vendorRoot = luminosa_miet_vendor_root();
            metalsFile = luminosa_miet_vendor_file('metals.mat');
            if isempty(vendorRoot)
                error('MIET-GUI was not found. Put it in external/miet-gui or vendor/miet-gui.');
            end
            if isempty(metalsFile)
                error('MIET-GUI support file metals.mat was not found in the vendored MIET tree.');
            end
            app.Metals = load(metalsFile);
            app.SpectrumFile = [];
            app.Monochrome = true;
        % place all "canvases" first:
            app.Figure = figure('MenuBar','none','Units','pixels',...           
                'Position',[100,100,12*app.widthIncr+140,16*app.heightIncr+100],...
                'NumberTitle','off','Name','Parameters of the sample',...
                'CloseRequestFcn',@app.closeWindow);
            app.PanelDispersion = uipanel('Parent',app.Figure,...
                'Title','Polychromatism','FontSize',10,...
                'BackgroundColor','white','Units','pixels','Position',...
                [20,12.5*app.heightIncr,12*app.widthIncr,2.5*app.heightIncr]);
            app.PanelUpper = uipanel('Parent',app.Figure,...
                'Title','Layers above molecule','FontSize',10,...
                'BackgroundColor','white','Units','pixels','Position',...
                [20,9*app.heightIncr,12*app.widthIncr,3*app.heightIncr]);
            app.PanelMolecule = uipanel('Parent',app.Figure,...
                'Title','Molecule''s layer','FontSize',10,...
                'BackgroundColor','white','Units','pixels','Position',...
                [20,4.5*app.heightIncr,5.5*app.widthIncr,4*app.heightIncr]);
            app.PanelGeneral = uipanel('Parent',app.Figure,...
                'Title','General parameters','FontSize',10,...
                'BackgroundColor','white','Units','pixels','Position',...
                [20+6.5*app.widthIncr,4.5*app.heightIncr,5.5*app.widthIncr,4*app.heightIncr]);
            app.PanelLower = uipanel('Parent',app.Figure,...
                'Title','Layers below molecule','FontSize',10,...
                'BackgroundColor','white','Units','pixels','Position',...
                [20,app.heightIncr,12*app.widthIncr,3*app.heightIncr]);
            
        % add control objects within those "canvases":
            app.DispText = uicontrol(app.PanelDispersion,'Units','Pixels','Position',...
                [app.widthIncr,5,10*app.widthIncr,2.5*app.heightIncr-25],...
                'Style','text','String',[sprintf('You can use only one wavelength ')...
                sprintf('for the evaluation (monochromatic) or load the emission ')...
                sprintf('spectrum of the dye and use all wavelengths allowed by ')...
                sprintf('your filters (polychromatic). In monocromatic mode, choose ')...
                sprintf('refractive indices according to wavelength. In polychromatic ')...
                sprintf('mode, use standard values for materials with low dispersion, ')...
                sprintf('e.g. glass or water, and placeholder values for metals.')]);
            app.DispChoose = uicontrol(app.Figure,'Style','popupmenu',...
                'String','monochromatic|polychromatic','Position',...
                [12*app.widthIncr+40,14.5*app.heightIncr-10,80,0.5*app.heightIncr],...
                'Callback',@app.ChooseDispersionType);
            app.PolyInfotext = uicontrol(app.Figure,'Units','Pixels','Position',...
                [12*app.widthIncr+40,9*app.heightIncr,80,5*app.heightIncr],...
                'Style','text','String',[sprintf('Placeholder values:\n')...
                sprintf('10=silver\n20=gold\n30=platinum\n40=paladium\n')...
                sprintf('50=copper\n60=aluminium\n70=chromium\n80=titanium\n')...
                sprintf('90=tungsten\n100=nickel\n110=beryllium\n120=ito')],...
                'Visible','off');
            app.RefrIndexUpper = uitable(app.PanelUpper,'Position',...
                [app.widthIncr,5,4*app.widthIncr,3*app.heightIncr-50],...
                'RowName',[],'ColumnName',[],'Data',[1.52;0],...
                'ColumnEditable',true,'ColumnWidth',{120},...
                'CellEditCallback',@app.manageTable);
            app.RefrIndexLower = uitable(app.PanelLower,'Position',...
                [app.widthIncr,5,4*app.widthIncr,3*app.heightIncr-50],...
                'RowName',[],'ColumnName',[],'Data',[0.3257+2.5792i;1.52;0],...
                'ColumnEditable',true,'ColumnWidth',{120},...
                'CellEditCallback',@app.manageTable);
            app.ThicknessUpper = uitable(app.PanelUpper,'Position',...
                [7*app.widthIncr,5,4*app.widthIncr,3*app.heightIncr-50],...
                'RowName',[],'ColumnName',[],'Data',0,...
                'ColumnEditable',true,'ColumnWidth',{120},...
                'CellEditCallback',@app.manageTable);
            app.ThicknessLower = uitable(app.PanelLower,'Position',...
                [7*app.widthIncr,5,4*app.widthIncr,3*app.heightIncr-50],...
                'RowName',[],'ColumnName',[],'Data',[12;0],...
                'ColumnEditable',true,'ColumnWidth',{120},...
                'CellEditCallback',@app.manageTable);
            app.RefrIndexMolec = uicontrol(app.PanelMolecule,'Position',...
                [app.widthIncr,2.45*app.heightIncr,4*app.widthIncr,20],...
                'Style','edit','String','1.52');
            app.ThicknessMolec = uicontrol(app.PanelMolecule,'Position',...
                [app.widthIncr,1.6*app.heightIncr,4*app.widthIncr,20],...
                'Style','edit','String','1000');
            app.QuantumYield = uicontrol(app.PanelMolecule,'Position',...
                [app.widthIncr,0.75*app.heightIncr,4*app.widthIncr,20],...
                'Style','edit','String','0.95');
            app.DipoleOrientation = uicontrol(app.PanelMolecule,'Position',...
                [app.widthIncr,0.05*app.heightIncr,4*app.widthIncr,20],...
                'Style','popupmenu','String','Fast rotating|Random fixed|Parallel|Vertical',...
                'Value',1);
            app.LaserMode = uicontrol(app.PanelGeneral,'Position',...
                [0.5*app.widthIncr,5*LayoutMolec+80,4*app.widthIncr,20],...
                'Style','popupmenu','String','radial|azimuthal|linear');
            app.NA = uicontrol(app.PanelGeneral,'Position',...
                [3.5*app.widthIncr,4*LayoutMolec+60,app.widthIncr,20],...
                'Style','edit','String','1.49');
            app.Focus = uicontrol(app.PanelGeneral,'Position',...
                [3.5*app.widthIncr,3*LayoutMolec+40,app.widthIncr,20],...
                'Style','edit','String','0');
            app.ExcWavelength = uicontrol(app.PanelGeneral,'Position',...
                [3.5*app.widthIncr,2*LayoutMolec+20,app.widthIncr,20],...
                'Style','edit','String','523');
            app.EmiWavelength = uicontrol(app.PanelGeneral,'Position',...
                [3.5*app.widthIncr,LayoutMolec,app.widthIncr,20],...
                'Style','edit','String','523','Callback',@app.FindRefrIndCallback);
            app.Wavel_Small = uicontrol(app.PanelGeneral,'Position',...
                [0.5*app.widthIncr,2*LayoutMolec+20,1.2*app.widthIncr,20],...
                'Style','edit','String','500','Visible','off');
            app.Wavel_Large = uicontrol(app.PanelGeneral,'Position',...
                [3.3*app.widthIncr,2*LayoutMolec+20,1.2*app.widthIncr,20],...
                'Style','edit','String','700','Visible','off');
            app.LoadSpectrumBt = uicontrol(app.Figure,'Style','pushbutton',...
                'String','Choose file','Units','pixels','Position',...
                [12*app.widthIncr+40,4.5*app.heightIncr+LayoutMolec,80,0.5*app.heightIncr],...
                'Visible','off','Callback',@app.ChooseFileCallback);
            app.SaveSettingsBt = uicontrol(app.Figure,'Style','pushbutton',...
                'String','Save Settings','Units','pixels',...
                'Position',[20,5,80,0.5*app.heightIncr],...
                'BusyAction','cancel','Callback',@app.SaveSettingsCallback);
            app.LoadSettingsBt = uicontrol(app.Figure,'Style','pushbutton',...
                'String','Load Settings','Units','pixels',...
                'Position',[110,5,80,0.5*app.heightIncr],...
                'BusyAction','cancel','Callback',@app.LoadSettingsCallback);
            app.ChooseMetal = uicontrol(app.Figure,'Style','popupmenu',...
                'String',['Silver|Gold|Platinum|Paladium|Copper|Aluminium'...
                '|Chromium|Titanium|Tungsten|Nickel|Beryllium|Ito'],...
                'Position',[290,5,70,0.5*app.heightIncr],...
                'Callback',@app.FindRefrIndCallback);
            app.ReadRefrInd = uicontrol(app.Figure,'Style','edit',...
                'Position',[370,5,100,0.5*app.heightIncr],'String','');
            app.OKButton = uicontrol(app.Figure,'Style','pushbutton',...
                'String','OK','Units','pixels',...
                'Position',[12*app.widthIncr+40,10+0.5*app.heightIncr,80,0.5*app.heightIncr],...
                'Callback',@app.OKButtonCallback);
            app.CancelButton = uicontrol(app.Figure,'Style','pushbutton',...
                'String','Cancel','Units','pixels',...
                'Position',[12*app.widthIncr+40,5,80,0.5*app.heightIncr],...
                'Callback',@app.CancelButtonCallback);
            app.ChooseCurveType = uicontrol(app.Figure,'Style','popupmenu',...
                'String','1st maximum|1st minimum','Position',...
                [12*app.widthIncr+40,6*app.heightIncr,80,0.5*app.heightIncr]);
            
        % add labels to the input fields
            app.Label_n1 = uicontrol(app.PanelUpper,'Position',...
                [app.widthIncr,3*app.heightIncr-40,4*app.widthIncr,18],...
                'Style','text','String','indices of refraction');
            app.Label_n1T = uicontrol(app.PanelUpper,'Position',...
                [5*app.widthIncr,3*app.heightIncr-68,2*app.widthIncr,18],...
                'Style','text','String','topmost');
            app.Label_n1B = uicontrol(app.PanelUpper,'Position',...
                [5*app.widthIncr,5,2*app.widthIncr,18],...
                'Style','text','String','bottommost');
            app.Label_d1 = uicontrol(app.PanelUpper,'Position',...
                [7*app.widthIncr,3*app.heightIncr-40,4*app.widthIncr,18],...
                'Style','text','String','thicknesses of layers [nm]');
            app.Label_n = uicontrol(app.PanelMolecule,'Position',...
                [app.widthIncr,3.05*app.heightIncr,4*app.widthIncr,20],...
                'Style','text','String','index of refraction');
            app.Label_d = uicontrol(app.PanelMolecule,'Position',...
                [app.widthIncr,2.2*app.heightIncr,4*app.widthIncr,20],...
                'Style','text','String','thickness of the layer [nm]');
            app.Label_QY = uicontrol(app.PanelMolecule,'Position',...
                [app.widthIncr,1.35*app.heightIncr,4*app.widthIncr,20],...
                'Style','text','String','free space quantum yield [0-1]');
            app.Label_DipoleOrientation = uicontrol(app.PanelMolecule,'Position',...
                [app.widthIncr,0.45*app.heightIncr,4*app.widthIncr,20],...
                'Style','text','String','dipole orientation');
            app.Label_n0 = uicontrol(app.PanelLower,'Position',...
                [app.widthIncr,3*app.heightIncr-40,4*app.widthIncr,18],...
                'Style','text','String','indices of refraction');
            app.Label_n0T = uicontrol(app.PanelLower,'Position',...
                [5*app.widthIncr,3*app.heightIncr-68,2*app.widthIncr,18],...
                'Style','text','String','topmost');
            app.Label_n0B = uicontrol(app.PanelLower,'Position',...
                [5*app.widthIncr,5,2*app.widthIncr,18],...
                'Style','text','String','bottommost');
            app.Label_d0 = uicontrol(app.PanelLower,'Position',...
                [7*app.widthIncr,4*app.heightIncr-40,4*app.widthIncr,18],...
                'Style','text','String','thicknesses of layers [nm]');
            app.Label_laserMode = uicontrol(app.PanelGeneral,'Style','text',...
                'Position',[0.5*app.widthIncr,6*LayoutMolec+100,4*app.widthIncr,40],...
                'String','polarization mode of laser');
            app.Label_NA = uicontrol(app.PanelGeneral,'Style','text',...
                'Position',[0.5*app.widthIncr,4*LayoutMolec+60,3*app.widthIncr,40],...
                'String','numerical aperture');
            app.Label_Focus = uicontrol(app.PanelGeneral,'Style','text',...
                'Position',[0.5*app.widthIncr,3*LayoutMolec+40,3*app.widthIncr,20],...
                'String','defocussing [nm]');
            app.Label_ExcWavelen = uicontrol(app.PanelGeneral,'Style','text',...
                'Position',[0.5*app.widthIncr,2*LayoutMolec+20,3*app.widthIncr,20],...
                'String','excitation wavel. [nm]');
            app.Label_EmiWavelen = uicontrol(app.PanelGeneral,'Style','text',...
                'Position',[0.5*app.widthIncr,LayoutMolec,3*app.widthIncr,20],...
                'String','emission wavel. [nm]');
            app.Label_WavelDisp = uicontrol(app.PanelGeneral,'Style','text',...
                'Position',[0.5*app.widthIncr,3*LayoutMolec+40,4*app.widthIncr,20],...
                'String','allowed wavelengths [nm]','Visible','off');
            app.Label_WavelRange = uicontrol(app.PanelGeneral,'Style','text',...
                'Position',[1.7*app.widthIncr,2*LayoutMolec+20,1.6*app.widthIncr,20],...
                'String','<=lambda<=','Visible','off');
            app.Label_SpectrumFile = uicontrol(app.PanelGeneral,'Style','text',...
                'Position',[0.5*app.widthIncr,LayoutMolec,4*app.widthIncr,20],...
                'String','File containing the spectrum...','Visible','off');
            app.Label_FindRefrInd = uicontrol(app.Figure,'Style','text',...
                'String','Find refr. index:','Units','pixels',...
                'Position',[200,5,80,0.5*app.heightIncr]);
            app.Label_ChooseCurveType = uicontrol(app.Figure,'Style','text',...
                'String','Calculate MIET calibration curve up to...','Position',...
                [12*app.widthIncr+40,6.7*app.heightIncr,80,1.5*app.heightIncr]);
            
            if nargin>=1 % if parameters have already been set in this session: use as new default
                input = varargin{1};
                if isfield(input,'NA') % was the correct structure passed as input argument?
                    set(app.RefrIndexUpper,'data',[input.RefrIndexUpper(end:-1:1) 0].');
                    set(app.RefrIndexLower,'data',[input.RefrIndexLower(end:-1:1) 0].');
                    set(app.ThicknessUpper,'data',[input.ThicknessUpper(end:-1:1) 0].');
                    set(app.ThicknessLower,'data',[input.ThicknessLower(end:-1:1) 0].');
                    set(app.RefrIndexMolec,'string',num2str(input.RefrIndexMolec));
                    set(app.ThicknessMolec,'string',num2str(input.ThicknessMolec));
                    set(app.QuantumYield,'string',num2str(input.QuantumYield));
                    set(app.LaserMode,'Value',input.LaserMode);
                    set(app.NA,'string',num2str(input.NA));
                    set(app.Focus,'string',num2str(input.Focus));
                    set(app.ExcWavelength,'string',num2str(input.ExcWavelength));
                    set(app.EmiWavelength,'string',num2str(input.EmiWavelength));
                    set(app.Wavel_Small,'string',num2str(input.Wavel_Small));
                    set(app.Wavel_Large,'string',num2str(input.Wavel_Large));
                    if ~isempty(input.SpectrumFile)
                        set(app.Label_SpectrumFile,'string',input.SpectrumFile);
                    end
                    app.SpectrumFile = input.SpectrumFile;
                    app.Monochrome = input.Monochrome;
                    if isfield(input,'DipoleOrientation')
                        set(app.DipoleOrientation,'Value',luminosa_miet_orientation_to_value(input.DipoleOrientation));
                    end
                    if app.Monochrome
                        set(app.DispChoose,'Value',1);
                    else
                        set(app.DispChoose,'Value',2);
                    end
                    if strcmp(input.CurveType,'maximum')
                        set(app.ChooseCurveType,'Value',1);
                    else
                        set(app.ChooseCurveType,'Value',2);
                    end
                end
            end
            ChooseDispersionType(app);
            FindRefrIndCallback(app); % update value in app.ReadRefrInd according to wavelength
            set(app.Figure,'ResizeFcn', @app.resizeWindow);
            app.resizeWindow(app.Figure, []);
        end

        function resizeWindow(app,src,~) % controls layout when window size changes
        % auxiliary constants
            fig = [];
            if ~isempty(app.Figure) && isgraphics(app.Figure)
                fig = app.Figure;
            elseif nargin >= 2 && ~isempty(src) && isgraphics(src)
                fig = src;
            else
                return;
            end
            figureSize = get(fig,'OuterPosition');
            if isempty(figureSize) || numel(figureSize) < 4
                figureSize = get(fig,'Position');
            end
            if isempty(figureSize) || numel(figureSize) < 4
                return;
            end
            app.widthIncr = (figureSize(3)-140)/12;
            app.heightIncr = figureSize(4)/16;
            LayoutMolec = (3.6*app.heightIncr-120)/7;
        % resize & redistribute "canvases"
            set(app.PanelDispersion,'Units','pixels','Position',...
                [20 12.5*app.heightIncr 12*app.widthIncr 2.5*app.heightIncr]);
            set(app.PanelUpper,'Units','pixels','Position',...
                [20 9*app.heightIncr 12*app.widthIncr 3*app.heightIncr]);
            set(app.PanelMolecule,'Units','pixels','Position',...
                [20 4.5*app.heightIncr 5.5*app.widthIncr 4*app.heightIncr]);
            set(app.PanelGeneral,'Units','pixels','Position',...
                [20+6.5*app.widthIncr,4.5*app.heightIncr,5.5*app.widthIncr,4*app.heightIncr]);
            set(app.PanelLower,'Units','pixels','Position',...
                [20,app.heightIncr,12*app.widthIncr,3*app.heightIncr]);
        % resize & redistribute control objects
            set(app.DispText,'Position',[app.widthIncr,5,10*app.widthIncr,2.5*app.heightIncr-25]);
            set(app.DispChoose,'Position',[12*app.widthIncr+40,14.5*app.heightIncr-10,80,0.5*app.heightIncr]);
            set(app.PolyInfotext,'Position',[12*app.widthIncr+40,9*app.heightIncr,80,5*app.heightIncr]);
            set(app.RefrIndexUpper,'Position',[app.widthIncr,5,4*app.widthIncr,3*app.heightIncr-50]);
            set(app.RefrIndexLower,'Position',[app.widthIncr,5,4*app.widthIncr,3*app.heightIncr-50]);
            set(app.ThicknessUpper,'Position',[7*app.widthIncr,5,4*app.widthIncr,3*app.heightIncr-50]);
            set(app.ThicknessLower,'Position',[7*app.widthIncr,5,4*app.widthIncr,3*app.heightIncr-50]);
            set(app.RefrIndexMolec,'Position',[app.widthIncr,2.45*app.heightIncr,4*app.widthIncr,20]);
            set(app.ThicknessMolec,'Position',[app.widthIncr,1.6*app.heightIncr,4*app.widthIncr,20]);
            set(app.QuantumYield,'Position',[app.widthIncr,0.75*app.heightIncr,4*app.widthIncr,20]);
            set(app.DipoleOrientation,'Position',[app.widthIncr,0.05*app.heightIncr,4*app.widthIncr,20]);
            set(app.EmiWavelength,'Position',[3.5*app.widthIncr,LayoutMolec,app.widthIncr,20]);
            set(app.Wavel_Small,'Position',[0.5*app.widthIncr,2*LayoutMolec+20,1.2*app.widthIncr,20]);
            set(app.Wavel_Large,'Position',[3.3*app.widthIncr,2*LayoutMolec+20,1.2*app.widthIncr,20]);
            set(app.LoadSpectrumBt,'Position',[12*app.widthIncr+40,4.5*app.heightIncr+LayoutMolec,80,0.5*app.heightIncr]);
            set(app.SaveSettingsBt,'Position',[20,5,80,0.5*app.heightIncr]);
            set(app.LoadSettingsBt,'Position',[110,5,80,0.5*app.heightIncr]);
            set(app.ChooseMetal,'Position',[290,5,70,0.5*app.heightIncr]);
            set(app.ReadRefrInd,'Position',[370,5,100,0.5*app.heightIncr]);
            set(app.OKButton,'Position',[12*app.widthIncr+40,10+0.5*app.heightIncr,80,0.5*app.heightIncr]);
            set(app.CancelButton,'Position',[12*app.widthIncr+40,5,80,0.5*app.heightIncr]);
            set(app.ChooseCurveType,'Position',[12*app.widthIncr+40,6*app.heightIncr,80,0.5*app.heightIncr]);
        % redistribute labels
            set(app.Label_n1,'Position',[app.widthIncr,3*app.heightIncr-40,4*app.widthIncr,18]);
            set(app.Label_n1T,'Position',[5*app.widthIncr,3*app.heightIncr-68,2*app.widthIncr,18]);
            set(app.Label_n1B,'Position',[5*app.widthIncr,5,2*app.widthIncr,18]);
            set(app.Label_d1,'Position',[7*app.widthIncr,3*app.heightIncr-40,4*app.widthIncr,18]);
            set(app.Label_n,'Position',[app.widthIncr,3.05*app.heightIncr,4*app.widthIncr,20]);
            set(app.Label_d,'Position',[app.widthIncr,2.2*app.heightIncr,4*app.widthIncr,20]);
            set(app.Label_QY,'Position',[app.widthIncr,1.35*app.heightIncr,4*app.widthIncr,20]);
            set(app.Label_DipoleOrientation,'Position',[app.widthIncr,0.45*app.heightIncr,4*app.widthIncr,20]);
            set(app.Label_n0,'Position',[app.widthIncr,3*app.heightIncr-40,4*app.widthIncr,18]);
            set(app.Label_n0T,'Position',[5*app.widthIncr,3*app.heightIncr-68,2*app.widthIncr,18]);
            set(app.Label_n0B,'Position',[5*app.widthIncr,5,2*app.widthIncr,18]);
            set(app.Label_d0,'Position',[7*app.widthIncr,3*app.heightIncr-40,4*app.widthIncr,18]);
            set(app.Label_EmiWavelen,'Position',[0.5*app.widthIncr,LayoutMolec,3*app.widthIncr,20]);
            set(app.Label_WavelDisp,'Position',[0.5*app.widthIncr,3*LayoutMolec+40,4*app.widthIncr,20]);
            set(app.Label_WavelRange,'Position',[1.7*app.widthIncr,2*LayoutMolec+20,1.6*app.widthIncr,20]);
            set(app.Label_SpectrumFile,'Position',[0.5*app.widthIncr,LayoutMolec,4*app.widthIncr,20]);
            set(app.Label_FindRefrInd,'Position',[200,5,80,0.5*app.heightIncr]);
            set(app.Label_ChooseCurveType,'Position',[12*app.widthIncr+40,6.7*app.heightIncr,80,1.5*app.heightIncr]);
        % things that are different for mono- and polychromatic mode
            disptype = get(app.DispChoose,'Value');
            switch disptype
                case 1 % monochromatic
                    set(app.LaserMode,'Position',[0.5*app.widthIncr,5*LayoutMolec+80,4*app.widthIncr,20]);
                    set(app.NA,'Position',[3.5*app.widthIncr,4*LayoutMolec+60,app.widthIncr,20]);
                    set(app.Focus,'Position',[3.5*app.widthIncr,3*LayoutMolec+40,app.widthIncr,20]);
                    set(app.ExcWavelength,'Position',[3.5*app.widthIncr,2*LayoutMolec+20,app.widthIncr,20]);
                    set(app.Label_laserMode,'Position',[0.5*app.widthIncr,6*LayoutMolec+100,4*app.widthIncr,20]);
                    set(app.Label_NA,'Position',[0.5*app.widthIncr,4*LayoutMolec+60,3*app.widthIncr,20]);
                    set(app.Label_Focus,'Position',[0.5*app.widthIncr,3*LayoutMolec+40,3*app.widthIncr,20]);
                    set(app.Label_ExcWavelen,'Position',[0.5*app.widthIncr,2*LayoutMolec+20,3*app.widthIncr,20]);
                case 2 % polychromatic
                    set(app.LaserMode,'Position',[2.5*app.widthIncr,6*LayoutMolec+100,2*app.widthIncr,20]);
                    set(app.NA,'Position',[1.3*app.widthIncr,5*LayoutMolec+80,app.widthIncr,20]);
                    set(app.Focus,'Position',[3.7*app.widthIncr,5*LayoutMolec+80,0.8*app.widthIncr,20]);
                    set(app.ExcWavelength,'Position',[3.5*app.widthIncr,4*LayoutMolec+60,app.widthIncr,20]);
                    set(app.Label_laserMode,'Position',[0.5*app.widthIncr,6*LayoutMolec+100,2*app.widthIncr,20]);
                    set(app.Label_NA,'Position',[0.5*app.widthIncr,5*LayoutMolec+80,0.8*app.widthIncr,20]);
                    set(app.Label_Focus,'Position',[2.3*app.widthIncr,5*LayoutMolec+80,1.4*app.widthIncr,20]);
                    set(app.Label_ExcWavelen,'Position',[0.5*app.widthIncr,4*LayoutMolec+60,3*app.widthIncr,20]);
            end
        end
      
        function manageTable(~,objectHandle,eventdata) % dynamically changes size of data entry matrices
            tabData = get(objectHandle,'data');
            if isnan(eventdata.NewData) % happens when you delete content of a cell
                temp = [tabData(1:eventdata.Indices(1)-1,:); tabData(eventdata.Indices(1)+1:end,:)];
                tabData = temp;
            end
            if tabData(end) ~= 0 % add a row filled with zeros at the end of the table
                tabData(end+1,:) = 0;
            end
            set(objectHandle,'data',tabData);
        end
          
        function SaveSettingsCallback(app,~,~) % save refr. indices, thicknesses, z-parameters
            [file,path] = uiputfile({'*.mat','MAT-files (*.mat)'},...
               'Please choose a filename for saving your settings');
            if file~=0
            % read setup parameters from the GUI
                monochrome=get(app.DispChoose,'Value'); % 1=monochromatic, 2=polychromatic
                n1=get(app.RefrIndexUpper,'data');  % refr. indices above molecule
                d1=get(app.ThicknessUpper,'data');  % thicknesses above molecule
                n0=get(app.RefrIndexLower,'data');  % refr. indices below molecule
                d0=get(app.ThicknessLower,'data');  % thicknesses below molecule
                n=get(app.RefrIndexMolec,'String'); % refr. index of molecule's layer
                d=get(app.ThicknessMolec,'String'); % thickness of molecule's layer
                QY=get(app.QuantumYield,'String');  % quantum yield of the molecule
                dipoleOrientation=luminosa_miet_orientation_from_value(get(app.DipoleOrientation,'Value'));
                laserMode=get(app.LaserMode,'Value'); % polarization mode of the excitation laser
                na=get(app.NA,'String');            % numerical aperture of objective
                focus=get(app.Focus,'String');      % from z_start to z_stop
                lambdaExc=get(app.ExcWavelength,'String');% excitation wavelength
                lambdaEmi=get(app.EmiWavelength,'String');% emission wavelength
                lambdaSmall=get(app.Wavel_Small,'String'); % smallest wavelength allowed by filters (polychromatic mode)
                lambdaLarge=get(app.Wavel_Large,'String'); % largest wavelength allowed by filters (polychromatic mode)
                curveType=get(app.ChooseCurveType,'Value');% calculate MIET calibration curve up to 1st maximum or last unique value?
                specFile=app.SpectrumFile;                 % file containing emission spectrum of emitter
                save([path file],'monochrome','n1','d1','n0','d0','n','d','QY','dipoleOrientation','laserMode','na',...
                    'focus','lambdaExc','lambdaEmi','lambdaSmall','lambdaLarge','curveType','specFile');
                disp([path file]);
            end
        end
        
        function LoadSettingsCallback(app,~,~) % load refr. indices, thicknesses, z-parameters
            [file,path] = uigetfile({'*.mat','MAT-files (*.mat)'},...
                'Please choose a filename for loading your settings.');
            if file~=0
                load([path file]);
                try
                    set(app.DispChoose,'Value',monochrome);
                    set(app.RefrIndexUpper,'data',n1);
                    set(app.ThicknessUpper,'data',d1);
                    set(app.RefrIndexLower,'data',n0);
                    set(app.ThicknessLower,'data',d0);
                    set(app.RefrIndexMolec,'String',n);
                    set(app.ThicknessMolec,'String',d);
                    set(app.QuantumYield,'String',QY);
                    if exist('dipoleOrientation','var')
                        set(app.DipoleOrientation,'Value',luminosa_miet_orientation_to_value(dipoleOrientation));
                    else
                        set(app.DipoleOrientation,'Value',1);
                    end
                    set(app.LaserMode,'Value',laserMode);
                    set(app.NA,'String',na);
                    set(app.Focus,'String',focus);
                    set(app.ExcWavelength,'String',lambdaExc);
                    set(app.EmiWavelength,'String',lambdaEmi);
                    set(app.Wavel_Small,'String',lambdaSmall);
                    set(app.Wavel_Large,'String',lambdaLarge);
                    set(app.ChooseCurveType,'Value',curveType);
                    set(app.Label_SpectrumFile,'String',specFile);
                    app.SpectrumFile = specFile;
                    FindRefrIndCallback(app);  % update value in app.ReadRefrInd according to wavelength
                    ChooseDispersionType(app); % update layout according to evaluation type (mono-/polychromatic)
                catch
                    warndlg('The specified file does not contain the correct kind of parameters.');
                end
            end
        end
       
        function ChooseDispersionType(app,~,~) % decide if mono- or polychromatic evaluation is used
            disptype = get(app.DispChoose,'Value');
            LayoutMolec = (3.6*app.heightIncr-120)/7;
            switch disptype
                case 1  % monochromatic
                    app.Monochrome = true;
                    
                    set(app.LaserMode,'Position',[0.5*app.widthIncr,5*LayoutMolec+80,4*app.widthIncr,20]);
                    set(app.NA,'Position',[3.5*app.widthIncr,4*LayoutMolec+60,app.widthIncr,20]);
                    set(app.Focus,'Position',[3.5*app.widthIncr,3*LayoutMolec+40,app.widthIncr,20]);
                    set(app.ExcWavelength,'Position',[3.5*app.widthIncr,2*LayoutMolec+20,app.widthIncr,20]);
                    set(app.Label_laserMode,'String','polarization mode of laser',...
                        'Position',[0.5*app.widthIncr,6*LayoutMolec+100,4*app.widthIncr,20]);
                    set(app.Label_NA,'String','numerical aperture','Position',[0.5*app.widthIncr,4*LayoutMolec+60,3*app.widthIncr,20]);
                    set(app.Label_Focus,'String','defocussing [nm]','Position',[0.5*app.widthIncr,3*LayoutMolec+40,3*app.widthIncr,20]);
                    set(app.Label_ExcWavelen,'Position',[0.5*app.widthIncr,2*LayoutMolec+20,3*app.widthIncr,20]);
                    
                    set(app.PolyInfotext,'Visible','off');
                    set(app.Label_WavelDisp,'Visible','off');
                    set(app.Wavel_Small,'Visible','off');
                    set(app.Label_WavelRange,'Visible','off');
                    set(app.Wavel_Large,'Visible','off');
                    set(app.LoadSpectrumBt,'Visible','off');
                    set(app.Label_SpectrumFile,'Visible','off');
                    set(app.EmiWavelength,'Visible','on');
                    set(app.Label_EmiWavelen,'Visible','on');
                    set(app.ChooseMetal,'Visible','on');
                    set(app.ReadRefrInd,'Visible','on');
                    set(app.Label_FindRefrInd,'Visible','on');
                    
                case 2 % polychromatic
                    app.Monochrome = false;

                    set(app.LaserMode,'Position',[2.5*app.widthIncr,6*LayoutMolec+100,2*app.widthIncr,20]);
                    set(app.NA,'Position',[1.3*app.widthIncr,5*LayoutMolec+80,app.widthIncr,20]);
                    set(app.Focus,'Position',[3.7*app.widthIncr,5*LayoutMolec+80,0.8*app.widthIncr,20]);
                    set(app.ExcWavelength,'Position',[3.5*app.widthIncr,4*LayoutMolec+60,app.widthIncr,20]);                    
                    set(app.Label_laserMode,'String','polarization',...
                        'Position',[0.5*app.widthIncr,6*LayoutMolec+100,2*app.widthIncr,20]);
                    set(app.Label_NA,'String','NA','Position',[0.5*app.widthIncr,5*LayoutMolec+80,0.8*app.widthIncr,20]);
                    set(app.Label_Focus,'String','defoc [nm]','Position',[2.3*app.widthIncr,5*LayoutMolec+80,1.4*app.widthIncr,20]);
                    set(app.Label_ExcWavelen,'Position',[0.5*app.widthIncr,4*LayoutMolec+60,3*app.widthIncr,20]);
                    
                    set(app.PolyInfotext,'Visible','on');
                    set(app.Label_WavelDisp,'Visible','on');
                    set(app.Wavel_Small,'Visible','on');
                    set(app.Label_WavelRange,'Visible','on');
                    set(app.Wavel_Large,'Visible','on');
                    set(app.LoadSpectrumBt,'Visible','on');
                    set(app.Label_SpectrumFile,'Visible','on');
                    set(app.EmiWavelength,'Visible','off');
                    set(app.Label_EmiWavelen,'Visible','off');
                    set(app.ChooseMetal,'Visible','off');
                    set(app.ReadRefrInd,'Visible','off');
                    set(app.Label_FindRefrInd,'Visible','off');
            end
        end
        
        function ChooseFileCallback(app,~,~) % choose the file containing the spectrum
            [file,path] = uigetfile({'*.asc','ASCII data (*.asc)';...
                '*.txt;*.dat','text files (*.txt, *.dat)';'*.*','All files (*.*)'},...
                'Please choose a file','MultiSelect','off');
            if file ~= 0
                app.SpectrumFile = [path file];
                set(app.Label_SpectrumFile,'String',app.SpectrumFile);
            end
        end
        
        function FindRefrIndCallback(app,~,~) % display the refr. index of a metal at a certain wavelength 
            lambda = str2double(get(app.EmiWavelength,'String'));
            metal = get(app.ChooseMetal,'Value');
            metalsFile = luminosa_miet_vendor_file('metals.mat');
            if isempty(metalsFile)
                set(app.ReadRefrInd,'String','');
                return;
            end
            metalsData = load(metalsFile, 'wavelength');
            wavelength = metalsData.wavelength;
            if 200<=lambda && lambda<=800 % wavelength range for which refr. ind. are saved
                switch metal % which metal was chosen by the user?
                    case 1
                        set(app.ReadRefrInd,'String',num2str(app.Metals.silver(wavelength==lambda)));
                    case 2
                        set(app.ReadRefrInd,'String',num2str(app.Metals.gold(wavelength==lambda)));
                    case 3
                        set(app.ReadRefrInd,'String',num2str(app.Metals.platinum(wavelength==lambda)));
                    case 4
                        set(app.ReadRefrInd,'String',num2str(app.Metals.palladium(wavelength==lambda)));
                    case 5
                        set(app.ReadRefrInd,'String',num2str(app.Metals.copper(wavelength==lambda)));
                    case 6
                        set(app.ReadRefrInd,'String',num2str(app.Metals.aluminum(wavelength==lambda)));
                    case 7
                        set(app.ReadRefrInd,'String',num2str(app.Metals.chromium(wavelength==lambda)));
                    case 8
                        set(app.ReadRefrInd,'String',num2str(app.Metals.titan(wavelength==lambda)));
                    case 9
                        set(app.ReadRefrInd,'String',num2str(app.Metals.tungsten(wavelength==lambda)));
                    case 10
                        set(app.ReadRefrInd,'String',num2str(app.Metals.nickel(wavelength==lambda)));
                    case 11
                        set(app.ReadRefrInd,'String',num2str(app.Metals.beryllium(wavelength==lambda)));
                    case 12
                        set(app.ReadRefrInd,'String',num2str(app.Metals.ito(wavelength==lambda)));
                end
            else
                set(app.ReadRefrInd,'String','no value available');
            end
        end
        
        function closeWindow(app,~,~) % runs when the window is closed
            delete(app.Figure)
        end
        
        function OKButtonCallback(~,~,~) % runs when button 'OK' is clicked
            uiresume;
        end
        
        function CancelButtonCallback(app,~,~) % runs when button 'cancel' is clicked
            delete(app.Figure)
        end
        
    end
end
