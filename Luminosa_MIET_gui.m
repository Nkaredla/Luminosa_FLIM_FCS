classdef Luminosa_MIET_gui < handle
% Graphical user interface to evaluate MIET data. Run this code to get an
% interface that will lead you through the evaluation.
%
% (c) Daja Ruhlandt, 2014
    
    properties
        Figure          % "canvas" for placing all other objects
        PanelInfo       % "canvas" for info-text
        PanelType       % "canvas" for choosing type of evaluation
        PanelFiles      % "canvas" for choosing files
        
        Infotext        % static text field containing information how to use the GUI
        EvaluationTypeRadioButtons % radio buttons to choose: MIET, smMIET, ...
        ChooseMIET      % --> MIET
        ChooseSmMIET    % --> smMIET
        MIETEvaluationType % popup menu to choose complete analysis or only show lifetime image
        smMIETEvaluationType % popup menu to choose complete analysis, only pattern matching, only show intensity image,...
        ParameterButton % pushbutton to enter sample parameters
        CutoffTime      % editable text field for entering cutoff time for calculation of average lifetime
        CutoffImage     % checkbox: use only photons arriving after cutoff for intensity image
        FreeLifetime    % editable text field for entering lifetime without presence of a cavity
        Button          % pushbutton to start main calculation
        ChooseFileButton % pushbutton to open dialog for choosing file (raw data or lifetimes)
        LTRadioButtons  % radiobuttons: do you want to calculate lifetime or use already calculated data?
        LTcalc          % -> calculate lifetimes
        LTuseOther      % -> use already calculated lifetimes
        FreeLTRadioButtons% radiobuttons: do you want to calculate free space lifetime or enter it directly?
        FreeLTcalc      % -> calculate free space lifetime
        FreeLTuseOther  % -> manually enter free space lifetime
        FreeLTButton    % pushbutton to choose file for calculating free space lifetime
        FreeLTCalcButton% pushbutton to calculate free space lifetime
        IRFRadioButtons % radiobuttons: do you have extra data for calculating IRF or do you want to estimate from normal data?
        IRFcalc         % -> calculate from extra data
        IRFestimate     % -> estimate from normal data
        IRFFileButton   % button to load file containing data for calculating IRF
        MakeImages      % checkbox: make nice images?
        MakeImageButton % pushbutton to plot already converted images
        Threshold       % editable text field for intensity threshold for images
        Binning         % editable text field for binning (bin width in px)
        
        Label_Lifetime  % label for lifetime without presence of a cavity
        Label_Threshold % label for intensity threshold for images
        Label_Binning   % label for binning
        Label_CutoffTime% label for cutoff time for calculation of average lifetime
        
        Label_Filename  % one of the .mat-files containing lifetimes or intensities
        Label_FreeSpaceLTFilename % .ht3 file containing image for calculation of free space lifetime
        Label_IRFFilename % .ht3 file containing data to calculate IRF
        
        Metals          % structure containing metal refractive indices
        Filename        % .mat-file containing already calculated lifetimes or intensities or rawdata
        IRFFilename     % .ht3-file containing data for calculating IRF
        FreeLTFilename  % .ht3-file containing data for calculating free space lifetime
        Parameters      % structure containing all sample parameters needed for evaluation
        Curve           % calibration curve
        Lifetime        % resulting lifetime image
        Height          % resulting height profile ('converted data')
        Intensity       % resulting intensity distribution
        heightIncr      % auxiliary variable for layout of app
        widthIncr       % auxiliary variable for layout of app
        LTcalcRaw       % auxiliary variable for remembering if LTcalc (true) or LTuseOther (fals)
    end

    
    methods
        % %  functions for designing and controlling the GUI  % % % % % % %
        function app = Luminosa_MIET_gui  % "constructor"
        % auxiliary constants (change if window size is changed, see resizeApp()) 
            app.heightIncr = 34;
            app.widthIncr = 34;
            smallHeightIncr = 0.34*app.heightIncr;
            app.LTcalcRaw = true;
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
            app.Filename='';
            app.IRFFilename=[];
            app.FreeLTFilename=[];
            app.Parameters=struct;
            app.Curve=[];
            app.Height=[];
            app.Intensity=[];
        % place all "canvases" first:
            app.Figure = figure('MenuBar','none','Units','pixels',...           
                'Position',[100,100,12*app.widthIncr+140+400,15*app.heightIncr+100],...
                'NumberTitle','off','Name','MIET GUI');
            app.PanelInfo = uipanel('Parent',app.Figure,...
                'Title','Info','FontSize',10,...
                'BackgroundColor','white','Units','pixels','Position',...
                [20,5.8*app.heightIncr+4*smallHeightIncr,12*app.widthIncr,4.8*app.heightIncr+6*smallHeightIncr]);
            app.PanelType = uipanel('Parent',app.Figure,...
                'Title','Type of evaluation','FontSize',10,...
                'BackgroundColor','white','Units','pixels','Position',...
                [20,app.heightIncr,12*app.widthIncr,3.6*app.heightIncr+5*smallHeightIncr]);
            app.PanelFiles = uipanel('Parent',app.Figure,'FontSize',10,...
                'Title','Choose source files and evaluation parameters',...
                'BackgroundColor','white','Units','pixels','Position',...
                [12*app.widthIncr+40,app.heightIncr,450,13*app.heightIncr]);
            
        % add control objects within those "canvases":
            app.Infotext = uicontrol(app.PanelInfo,'Units','Pixels','Position',...
                [app.widthIncr,smallHeightIncr,10*app.widthIncr,4.8*app.heightIncr+3*smallHeightIncr],...
                'Style','text','String',[sprintf('How to use this GUI:\n')...
                sprintf('1. Choose the evaluation mode: standard MIET (pixel-by-')...
                sprintf('pixel) or ')...
                sprintf('single molecule MIET (using pattern matching).\n')...
                sprintf('2. Press the button ''Define sample parameters'' to ')...
                sprintf('set parameters such as wavelength, refractive indices etc.\n')...
                sprintf('3. Choose the evaluation submode: only lifetime and ')...
                sprintf('intensity image, MIET height profile, different fitting ')...
                sprintf('modes etc.\n')...
                sprintf('4. Only for standard MIET: estimate how many ns after the peak ')...
                sprintf('of the tcspc-curve the curve looks like an exponential decay ')...
                sprintf('("cutoff-time").\n')...
                sprintf('5. Only for standard MIET: decide if you want to read in ')...
                sprintf('an unprocessed ht3-file or if you want to evaluate data ')...
                sprintf('that has already been converted to lifetimes (needs the ')...
                sprintf('file ending on _PS.mat).\n')...
                sprintf('6. Only for single-molecule MIET: decide if you ')...
                sprintf('want to use another file containing ')...
                sprintf('experimental data to find the instrument response function ')...
                sprintf('or if you want to estimate it from your normal data.\n')...
                sprintf('7. Use a third filename to calculate the free space ')...
                sprintf('lifetime of your dye via three different methods or ')...
                sprintf('directly type in the value if it is known to you.\n')...
                sprintf('8. If you like, you can visualize the height profile in 3D.')]);
            app.EvaluationTypeRadioButtons = uibuttongroup(app.PanelType,...
                'Units','pixels','Visible','on','Position',...
                [app.widthIncr,2.4*app.heightIncr+3*smallHeightIncr,10*app.widthIncr,1.2*app.heightIncr],...
                'SelectionChangeFcn',@app.EvaluationTypeCallback);
            app.ChooseMIET = uicontrol(app.EvaluationTypeRadioButtons,...
                'String','Evaluate MIET data pixel-by-pixel','Units','pixels',...
                'Position',[10,0.65*app.heightIncr,250,0.45*app.heightIncr],...
                'Style','radiobutton','Max',1,'Min',0,'Value',1);
            app.ChooseSmMIET = uicontrol(app.EvaluationTypeRadioButtons,...
                'String','Use single-molecule patterns to evaluate MIET data','Units',...
                'pixels','Style','radiobutton','Max',1,'Min',0,'Value',0,...
                'Position',[10,0.1*app.heightIncr,250,0.45*app.heightIncr]);
            app.MIETEvaluationType = uicontrol(app.PanelType,'Style','popupmenu',...
                'String',['Complete analysis|Only show lifetime and intensity images'...
                '|Only show MIET lifetime vs. height calibration curve'],...
                'Units','pixels','Visible','on','Position',...
                [app.widthIncr,1.05*app.heightIncr,10*app.widthIncr,0.5*app.heightIncr]);
            app.smMIETEvaluationType = uicontrol(app.PanelType,'Style','popupmenu',...
                'String',['Complete analysis of single molecules using pattern detection'...
                '|Complete analysis of elliptical regions of interest'...
                '|Only lifetime and intensity using pattern detection, no height'...
                '|Lifetime and intensity of randomly oriented molecules, e.g. fluorescent beads'...
                '|Lifetime and intensity of elliptical regions of interest'...
                '|Find and display patterns and raw intensity image'...
                '|Only show raw intensity image'...
                '|Only show MIET lifetime vs. height calibration curves'],...
                'Units','pixels','Visible','off','Position',...
                [app.widthIncr,1.05*app.heightIncr,10*app.widthIncr,0.5*app.heightIncr]);
            app.ParameterButton = uicontrol(app.PanelType,'Style','pushbutton',...
                'String','Define sample parameters','Units','pixels','Position',...
                [app.widthIncr,1.55*app.heightIncr+2*smallHeightIncr,10*app.widthIncr,0.5*app.heightIncr],...
                'Callback',@app.SampleParametersCallback,'BusyAction','cancel');
            app.CutoffTime = uicontrol(app.PanelType,'Style','edit',...
                'String','0.5','Units','pixels','Position',...
                [3.2*app.widthIncr,0.3*app.heightIncr+0.75*smallHeightIncr,0.8*app.widthIncr,0.5*app.heightIncr]);
            app.CutoffImage = uicontrol(app.PanelType,'Units','pixels','Position',...
                [4.4*app.widthIncr,0.3*app.heightIncr+0.75*smallHeightIncr,6.6*app.widthIncr,0.5*app.heightIncr],...
                'Style','checkbox','String','Only photons after cutoff in intensity image?',...
                'Min',0,'Max',1,'Value',0);
            app.Button = uicontrol(app.Figure,'Style','pushbutton',...
                'String','Evaluate','Units','pixels',...
                'Position',[12*app.widthIncr+40,5,80,0.5*app.heightIncr],...
                'Callback',@app.btnCallback,'BusyAction','cancel');
            app.ChooseFileButton = uicontrol(app.PanelFiles,'Style','pushbutton',...
                'String','Choose file','Units','pixels','Visible','on',...
                'Position',[320,5.7*app.heightIncr,100,0.5*app.heightIncr],...
                'Callback',@app.ChooseFileBtnCallback,'BusyAction','cancel');
            app.LTRadioButtons = uibuttongroup(app.PanelFiles,'Units','pixels',...
                'Position',[10,8.5*app.heightIncr-60,300,1.5*app.heightIncr],'Visible','on');
            app.LTcalc = uicontrol(app.LTRadioButtons,...
                'String','Calculate lifetimes from raw data [.ht3-file]',...
                'Position',[10,0.65*app.heightIncr,250,0.45*app.heightIncr],...
                'Style','radiobutton','Max',1,'Min',0,'Value',1);
            app.LTuseOther = uicontrol(app.LTRadioButtons,...
                'String','Use previously calculated lifetimes [.mat-file]',...
                'Position',[10,0.1*app.heightIncr,250,0.45*app.heightIncr],...
                'Style','radiobutton','Max',1,'Min',0,'Value',0);
            app.FreeLTRadioButtons = uibuttongroup(app.PanelFiles,'Units','pixels',...
                'Position',[10,6*app.heightIncr+6*smallHeightIncr,300,1.2*app.heightIncr],...
                'SelectionChangeFcn',@app.FreeSpaceLifetimeTypeCallback,'Visible','on');      
            app.FreeLTcalc = uicontrol(app.FreeLTRadioButtons,...
                'String','Calculate free space lifetime from raw data',...
                'Position',[10,0.1*app.heightIncr,250,0.45*app.heightIncr],...
                'Style','radiobutton','Max',1,'Min',0,'Value',0);
            app.FreeLTuseOther = uicontrol(app.FreeLTRadioButtons,...
                'String','Manually enter free space lifetime [ns]',...
                'Position',[10,0.65*app.heightIncr,250,0.45*app.heightIncr],...
                'Style','radiobutton','Max',1,'Min',0,'Value',1);
            app.FreeLifetime = uicontrol(app.PanelFiles,'Style','edit',...
                'Position',[320,3.45*app.heightIncr,100,0.5*app.heightIncr],...
                'String','4.3','Enable','on');
            app.FreeLTButton = uicontrol(app.PanelFiles,'Style','pushbutton',...
                'String','Choose file','Units','pixels','Visible','off',...
                'Position',[320,2.6*app.heightIncr,100,0.5*app.heightIncr],...
                'Callback',@app.ChooseFileBtnCallback,'BusyAction','cancel');
            app.FreeLTCalcButton = uicontrol(app.PanelFiles,'Style','pushbutton',...
                'String','Calculate LT','Units','pixels','Visible','off',...
                'Position',[320,1.9*app.heightIncr,100,0.5*app.heightIncr],...
                'Callback',@app.GetFreeSpaceLifetime,'BusyAction','cancel');
            app.IRFRadioButtons = uibuttongroup(app.PanelFiles,'Units','pixels',...
                'Position',[10,3.6*app.heightIncr+4*smallHeightIncr,300,1.2*app.heightIncr],...
                'Visible','off','SelectionChangeFcn',@app.IRFDeterminationTypeCallback);
            app.IRFcalc = uicontrol(app.IRFRadioButtons,...
                'String','Calculate IRF from extra data','Position',...
                [10,0.1*app.heightIncr,250,0.45*app.heightIncr],...
                'Style','radiobutton','Max',1,'Min',0,'Value',0);
            app.IRFestimate = uicontrol(app.IRFRadioButtons,...
                'String','Estimate IRF','Position',...
                [10,0.65*app.heightIncr,250,0.45*app.heightIncr],...
                'Style','radiobutton','Max',1,'Min',0,'Value',1);
            app.IRFFileButton = uicontrol(app.PanelFiles,'Style','pushbutton',...
                'String','Choose file','Units','pixels','Visible','off',...
                'Position',[320,2.4*app.heightIncr+3*smallHeightIncr,100,0.5*app.heightIncr],...
                'Callback',@app.ChooseFileBtnCallback,'BusyAction','cancel');
            app.MakeImages = uicontrol(app.PanelFiles,'Units','pixels',...
                'Position',[10,0.7*app.heightIncr+smallHeightIncr,300,0.5*app.heightIncr],...
                'Style','checkbox','String','Visualize height profile?',...
                'Min',0,'Max',1,'Value',0,'Callback',@app.ImageCheckboxCallback);
            app.MakeImageButton = uicontrol(app.PanelFiles,'Style','pushbutton',...
                'String','Visualize now','Units','pixels','Visible','off',...
                'Position',[320,2*app.heightIncr,100,0.5*app.heightIncr],...
                'Callback',@app.MakeImageBtnCallback,'BusyAction','cancel');
            app.Threshold = uicontrol(app.PanelFiles,'Units','pixels',...
                'Position',[320,app.heightIncr,100,20],'Style','edit',...
                'String','50','Visible','off');
            app.Binning = uicontrol(app.PanelType,'Style','edit',...
                'String','1','Visible','off','Position',...
                [9*app.widthIncr,0.7*app.heightIncr+smallHeightIncr,2*app.widthIncr,0.5*app.heightIncr]);
            
        % add labels to the input fields
            app.Label_Filename = uicontrol(app.PanelFiles,'Units','pixels',...
                'Style','text','String','Please choose the truncated name of the file...',...
                'Position',[10,5*app.heightIncr,300,40],'Visible','on');
            app.Label_Lifetime = uicontrol(app.PanelFiles,'Style','text',...
                'Position',[320,4.15*app.heightIncr,100,0.5*app.heightIncr],...
                'String','free space LT [ns]');
            app.Label_FreeSpaceLTFilename = uicontrol(app.PanelFiles,'Units','pixels',...
                'Style','text','String','Please choose the file used to calculate free space lifetime...',...
                'Position',[10,1.9*app.heightIncr,300,1.2*app.heightIncr],'Visible','off');
            app.Label_IRFFilename = uicontrol(app.PanelFiles,'Units','pixels',...
                'Style','text','String','Please choose the name of the file containing the IRF data...',...
                'Position',[10,2.4*app.heightIncr+3*smallHeightIncr,300,1.2*app.heightIncr],'Visible','off');
            app.Label_Threshold = uicontrol(app.PanelFiles,'Units','pixels',...
                'Style','text','String','Lower intensity threshold for images [counts]',...
                'Position',[10,0.35*app.heightIncr,300,0.5*app.heightIncr],'Visible','off');
            app.Label_Binning = uicontrol(app.PanelType,'Units','pixels',...
                'Style','text','String','Binning: bin width[px]','Visible','off',...
                'Position',[app.widthIncr,1.05*app.heightIncr,7.5*app.widthIncr,0.5*app.heightIncr]);
            app.Label_CutoffTime = uicontrol(app.PanelType,'Units','pixels',...
                'Style','text','String','Cutoff-time [ns]','Position',...
                [app.widthIncr,0.3*app.heightIncr+0.75*smallHeightIncr,2.2*app.widthIncr,0.5*app.heightIncr]);
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
            app.widthIncr = (figureSize(3)-140-400)/12;
            app.heightIncr = figureSize(4)/15;
            smallHeightIncr = 0.34*app.heightIncr;
        % resize & redistribute "canvases"
            set(app.PanelInfo,'Position',[20,5.8*app.heightIncr+4*smallHeightIncr,12*app.widthIncr,4.8*app.heightIncr+6*smallHeightIncr]);
            set(app.PanelType,'Position',[20,app.heightIncr,12*app.widthIncr,3.6*app.heightIncr+5*smallHeightIncr]);
            set(app.PanelFiles,'Units','pixels','Position',...
                [12*app.widthIncr+40,app.heightIncr,450,13*app.heightIncr]);
        % resize & redistribute control objects
            set(app.Infotext,'Position',[app.widthIncr,smallHeightIncr,10*app.widthIncr,4.8*app.heightIncr+3*smallHeightIncr]);
            set(app.Button,'Position',[12*app.widthIncr+40+370,5,80,0.5*app.heightIncr]);
            set(app.EvaluationTypeRadioButtons,'Position',[app.widthIncr,2.4*app.heightIncr+3*smallHeightIncr,10*app.widthIncr,1.2*app.heightIncr]);
            set(app.ChooseMIET,'Position',[0.5*app.widthIncr,0.65*app.heightIncr,8*app.widthIncr,0.45*app.heightIncr]);
            set(app.ChooseSmMIET,'Position',[0.5*app.widthIncr,0.1*app.heightIncr,8*app.widthIncr,0.45*app.heightIncr]);
            set(app.MIETEvaluationType,'Position',[app.widthIncr,app.heightIncr+1.5*smallHeightIncr,10*app.widthIncr,0.5*app.heightIncr]);
            set(app.smMIETEvaluationType,'Position',[app.widthIncr,app.heightIncr+1.5*smallHeightIncr,10*app.widthIncr,0.5*app.heightIncr]);
            set(app.ParameterButton,'Position',[app.widthIncr,1.7*app.heightIncr+2.25*smallHeightIncr,10*app.widthIncr,0.5*app.heightIncr]);
            set(app.CutoffTime,'Position',[3.2*app.widthIncr,0.3*app.heightIncr+0.75*smallHeightIncr,0.8*app.widthIncr,0.5*app.heightIncr]);
            set(app.CutoffImage,'Position',[4.4*app.widthIncr,0.3*app.heightIncr+0.75*smallHeightIncr,6.6*app.widthIncr,0.5*app.heightIncr]);
            set(app.ChooseFileButton,'Position',[320,7.2*app.heightIncr+7*smallHeightIncr,100,0.5*app.heightIncr]);
            set(app.LTRadioButtons,'Position',[10,8.4*app.heightIncr+8*smallHeightIncr,300,1.2*app.heightIncr]);
            set(app.LTcalc,'Position',[10,0.65*app.heightIncr,250,0.45*app.heightIncr]);
            set(app.LTuseOther,'Position',[10,0.1*app.heightIncr,250,0.45*app.heightIncr]);
            set(app.FreeLTcalc,'Position',[10,0.1*app.heightIncr,250,0.45*app.heightIncr]);
            set(app.FreeLTuseOther,'Position',[10,0.65*app.heightIncr,250,0.45*app.heightIncr]);
            set(app.FreeLTRadioButtons,'Position',[10,3.6*app.heightIncr+4*smallHeightIncr,300,1.2*app.heightIncr]);
            set(app.FreeLifetime,'Position',[320,3.6*app.heightIncr+4*smallHeightIncr,100,0.5*app.heightIncr]);
            set(app.FreeLTButton,'Position',[320,3.1*app.heightIncr+3*smallHeightIncr,100,0.5*app.heightIncr]);
            set(app.FreeLTCalcButton,'Position',[320,2.4*app.heightIncr+3*smallHeightIncr,100,0.5*app.heightIncr]);
            set(app.IRFRadioButtons,'Position',[10,6*app.heightIncr+6*smallHeightIncr,300,1.2*app.heightIncr]);
            set(app.IRFcalc,'Position',[10,0.1*app.heightIncr,250,0.45*app.heightIncr]);
            set(app.IRFestimate,'Position',[10,0.65*app.heightIncr,250,0.45*app.heightIncr]);
            set(app.IRFFileButton,'Position',[320,4.8*app.heightIncr+5*smallHeightIncr,100,0.5*app.heightIncr]);
            set(app.MakeImages,'Position',[10,0.7*app.heightIncr+smallHeightIncr,300,0.5*app.heightIncr]);
            set(app.MakeImageButton,'Position',[320,0.7*app.heightIncr+smallHeightIncr,100,0.5*app.heightIncr]);
            set(app.Threshold,'Position',[320,smallHeightIncr,100,0.5*app.heightIncr]);
            set(app.Binning,'Position',[9*app.widthIncr,0.3*app.heightIncr+0.75*smallHeightIncr,2*app.widthIncr,0.5*app.heightIncr]);
        % redistribute labels
            set(app.Label_Filename,'Position',[10,7.2*app.heightIncr+7*smallHeightIncr,300,1.2*app.heightIncr]);
            set(app.Label_Lifetime,'Position',[320,4.3*app.heightIncr+4*smallHeightIncr,100,0.5*app.heightIncr]);
            set(app.Label_FreeSpaceLTFilename,'Position',[10,2.4*app.heightIncr+3*smallHeightIncr,300,1.2*app.heightIncr]);
            set(app.Label_IRFFilename,'Position',[10,4.8*app.heightIncr+5*smallHeightIncr,300,1.2*app.heightIncr]);
            set(app.Label_Threshold,'Position',[10,smallHeightIncr,300,0.5*app.heightIncr]);
            set(app.Label_Binning,'Position',[app.widthIncr,0.3*app.heightIncr+0.75*smallHeightIncr,7.5*app.widthIncr,0.5*app.heightIncr]);
            set(app.Label_CutoffTime,'Position',[app.widthIncr,0.3*app.heightIncr+0.75*smallHeightIncr,2.2*app.widthIncr,0.5*app.heightIncr]);
        end
        
        function closeWindow(app,~,~) % runs when the window is closed
            delete(app.Figure)
        end
        
        
        % %  functions for choosing settings, files,... % % % % % % % % % %
        function EvaluationTypeCallback(app,~,eventdata)
            if eventdata.NewValue == app.ChooseMIET
                set(app.ParameterButton,'Callback',@app.SampleParametersCallback);
                set(app.Button,'Callback',@app.btnCallback);
                set(app.Binning,'Visible','off');
                set(app.Label_Binning,'Visible','off');
                set(app.CutoffTime,'Visible','on');
                set(app.Label_CutoffTime,'Visible','on');
                set(app.CutoffImage,'Visible','on');
                set(app.MIETEvaluationType,'Visible','on');
                set(app.smMIETEvaluationType,'Visible','off');
                set(app.LTuseOther,'Visible','on');
                if app.LTcalcRaw
                    set(app.LTcalc,'Value',1); set(app.LTuseOther,'Value',0);
                else
                    set(app.LTcalc,'Value',0); set(app.LTuseOther,'Value',1);
                end
                set(app.IRFRadioButtons,'Visible','off');
                set(app.IRFFileButton,'Visible','off');
                set(app.Label_IRFFilename,'Visible','off');
                set(app.MakeImages,'Visible','on');
                if get(app.MakeImages,'Value')==get(app.MakeImages,'Max')
                    set(app.MakeImageButton,'Visible','on');
                    set(app.Threshold,'Visible','on');
                    set(app.Label_Threshold,'Visible','on');
                end
            else
                app.LTcalcRaw = get(app.LTcalc,'Value')==1;
                set(app.ParameterButton,'Callback',@app.smMIETParametersCallback);
                set(app.Button,'Callback',@app.smMIETCallback);
                set(app.Binning,'Visible','off');
                set(app.Label_Binning,'Visible','off');
                set(app.CutoffTime,'Visible','off');
                set(app.Label_CutoffTime,'Visible','off');
                set(app.CutoffImage,'Visible','off');
                set(app.MIETEvaluationType,'Visible','off');
                set(app.smMIETEvaluationType,'Visible','on');
                set(app.LTuseOther,'Visible','off');
                set(app.LTcalc,'Value',1);
                set(app.IRFRadioButtons,'Visible','on');
                if(get(app.IRFcalc,'Value')==1)
                    set(app.IRFFileButton,'Visible','on');
                    set(app.Label_IRFFilename,'Visible','on');
                end
                set(app.MakeImages,'Visible','off');
                set(app.MakeImageButton,'Visible','off');
                set(app.Threshold,'Visible','off');
                set(app.Label_Threshold,'Visible','off');
            end
        end
        
        function SampleParametersCallback(app,~,~) % open new window for entering sample parameters
            if isfield(app.Parameters,'Wavelength')
                handle = Luminosa_MIET_parameters(app.Parameters);   % handle of the new window containing already set parameters
            else
                handle = Luminosa_MIET_parameters;   % handle of the new window with default parameters
            end
            uiwait(handle.Figure);      % block until new window is closed or until command 'uiresume'
        % if user clicked 'cancel' or closed window, do not change parameters;
        % if user clicked 'OK', all the data is accessible via handle and window still has to be closed   
            try
                tempdata=get(handle.RefrIndexUpper,'data'); % refr. indices above molecule
                app.Parameters.RefrIndexUpper = []; % reset array
                app.Parameters.RefrIndexUpper(numel(tempdata)-1:-1:1)=tempdata(1:end-1);
                tempdata=get(handle.RefrIndexLower,'data'); % refr. indices below molecule
                app.Parameters.RefrIndexLower = []; % reset array
                app.Parameters.RefrIndexLower(numel(tempdata)-1:-1:1)=tempdata(1:end-1);
                tempdata=get(handle.ThicknessUpper,'data'); % thicknesses above molecule
                app.Parameters.ThicknessUpper = []; % reset array
                app.Parameters.ThicknessUpper(numel(tempdata)-1:-1:1)=tempdata(1:end-1);
                tempdata=get(handle.ThicknessLower,'data'); % thicknesses below molecule
                app.Parameters.ThicknessLower = []; % reset array
                app.Parameters.ThicknessLower(numel(tempdata)-1:-1:1)=tempdata(1:end-1);
                app.Parameters.RefrIndexMolec = str2double(get(handle.RefrIndexMolec,'String'));
                app.Parameters.ThicknessMolec = str2double(get(handle.ThicknessMolec,'String'));
                app.Parameters.Wavelength     = str2double(get(handle.Wavelength,'String'));
                app.Parameters.Z_Start        = str2double(get(handle.Z_Start,'String'));
                app.Parameters.Z_Stop         = str2double(get(handle.Z_Stop,'String'));
                app.Parameters.Z_NumSteps     = str2double(get(handle.Z_NumSteps,'String'));
                app.Parameters.Wavel_Small    = str2double(get(handle.Wavel_Small,'String'));
                app.Parameters.Wavel_Large    = str2double(get(handle.Wavel_Large,'String'));
                app.Parameters.SpectrumFile   = handle.SpectrumFile;
                app.Parameters.Monochrome     = handle.Monochrome;
                app.Parameters.QuantumYield   = str2double(get(handle.QuantumYield,'String'));
                app.Parameters.DipoleOrientation = luminosa_miet_orientation_from_value(get(handle.DipoleOrientation,'Value'));
                switch get(handle.ChooseCurveType,'Value')
                    case 1; app.Parameters.CurveType = 'maximum';
                    case 2; app.Parameters.CurveType = 'minimum';
                    case 3; app.Parameters.CurveType = 'manual';
                end
                close 'Parameters of the sample'    % close the other window
                disp('Successfully changed sample parameters.');
            catch
                disp('Setting parameters aborted - nothing was changed.');
            end
        end
        
        function smMIETParametersCallback(app,~,~) % open new window for entering smMIET parameters
            if isfield(app.Parameters,'NA')
                handle = Luminosa_smMIET_parameters(app.Parameters);   % handle of the new window containing already set parameters
            else
                handle = Luminosa_smMIET_parameters;   % handle of the new window with default parameters
            end
            uiwait(handle.Figure);      % block until new window is closed or until command 'uiresume'
        % if user clicked 'cancel' or closed window, do not change parameters;
        % if user clicked 'OK', all the data is accessible via handle and window still has to be closed   
            try
                tempdata=get(handle.RefrIndexUpper,'data'); % refr. indices above molecule
                app.Parameters.RefrIndexUpper = []; % reset array
                app.Parameters.RefrIndexUpper(numel(tempdata)-1:-1:1)=tempdata(1:end-1);
                tempdata=get(handle.RefrIndexLower,'data'); % refr. indices below molecule
                app.Parameters.RefrIndexLower = []; % reset array
                app.Parameters.RefrIndexLower(numel(tempdata)-1:-1:1)=tempdata(1:end-1);
                tempdata=get(handle.ThicknessUpper,'data'); % thicknesses above molecule
                app.Parameters.ThicknessUpper = []; % reset array
                app.Parameters.ThicknessUpper(numel(tempdata)-1:-1:1)=tempdata(1:end-1);
                tempdata=get(handle.ThicknessLower,'data'); % thicknesses below molecule
                app.Parameters.ThicknessLower = []; % reset array
                app.Parameters.ThicknessLower(numel(tempdata)-1:-1:1)=tempdata(1:end-1);
                app.Parameters.RefrIndexMolec = str2double(get(handle.RefrIndexMolec,'String'));
                app.Parameters.ThicknessMolec = str2double(get(handle.ThicknessMolec,'String'));
                app.Parameters.LaserMode      = get(handle.LaserMode,'Value');
                app.Parameters.NA             = str2double(get(handle.NA,'String'));
                app.Parameters.Focus          = str2double(get(handle.Focus,'String'));
                app.Parameters.ExcWavelength  = str2double(get(handle.ExcWavelength,'String'));
                app.Parameters.EmiWavelength  = str2double(get(handle.EmiWavelength,'String'));
                app.Parameters.Wavel_Small    = str2double(get(handle.Wavel_Small,'String'));
                app.Parameters.Wavel_Large    = str2double(get(handle.Wavel_Large,'String'));
                app.Parameters.SpectrumFile   = handle.SpectrumFile;
                app.Parameters.Monochrome     = handle.Monochrome;
                app.Parameters.QuantumYield   = str2double(get(handle.QuantumYield,'String'));
                app.Parameters.DipoleOrientation = luminosa_miet_orientation_from_value(get(handle.DipoleOrientation,'Value'));
                switch get(handle.ChooseCurveType,'Value')
                    case 1; app.Parameters.CurveType = 'maximum';
                    case 2; app.Parameters.CurveType = 'minimum';
                end
                close 'Parameters of the sample'    % close the other window
                disp('Successfully changed sample parameters.');
            catch
                disp('Setting parameters aborted - nothing was changed.');
            end
        end
        
        function FreeSpaceLifetimeTypeCallback(app,~,eventdata) % controls free space LT radiobuttons
            if eventdata.NewValue == app.FreeLTcalc % evaluate raw data?
                set(app.Label_FreeSpaceLTFilename,'Visible','on');
                set(app.FreeLTButton,'Visible','on');
                set(app.FreeLTCalcButton,'Visible','on');
            else                    % manually type in free space lifetime
                set(app.Label_FreeSpaceLTFilename,'Visible','off');
                set(app.FreeLTButton,'Visible','off');
                set(app.FreeLTCalcButton,'Visible','off');
            end
        end
        
        function IRFDeterminationTypeCallback(app,~,eventdata) % controls IRF radiobuttons
            if eventdata.NewValue == app.IRFcalc % calculate IRF from extra data?
                set(app.Label_IRFFilename,'Visible','on');
                set(app.IRFFileButton,'Visible','on');
            else
                set(app.Label_IRFFilename,'Visible','off');
                set(app.IRFFileButton,'Visible','off');
            end
        end
        
        function ImageCheckboxCallback(app,objectHandle,~) % runs when visualization checkbox is clicked
            if (get(objectHandle,'Value') == get(objectHandle,'Max'))
                set(app.Label_Threshold,'Visible','on');
                set(app.Threshold,'Visible','on');
                set(app.MakeImageButton,'Visible','on');
            else
                set(app.Label_Threshold,'Visible','off');
                set(app.Threshold,'Visible','off');
                set(app.MakeImageButton,'Visible','off');
            end
        end
                
        function ChooseFileBtnCallback(app,objectHandle,~) % runs when one of the "chooseFile" buttons is clicked
            useRaw = get(app.LTcalc,'Value')==1;    % use raw data?
            LTfile = objectHandle==app.ChooseFileButton; % is this button for the lifetime?
            if LTfile && ~useRaw % use previously calculated lifetimes
                [file,path] = uigetfile({'*.mat','MAT-files (*.mat)';...
                    '*.*','All files (*.*)'},'Please choose a file','MultiSelect','on');
            else
                [file,path] = uigetfile({'*.ht3;*.ptu','Raw data (*.ht3,*.ptu)';...
                    '*.mat','MAT-files (*.mat)'; '*.*','All files (*.*)'},...
                    'Please choose a file','MultiSelect','on');
            end
            if iscell(file)
                filename=cell(size(file));
                for i=1:size(file,2)
                % raw data is in .ht3-file, user-defined data may have any name
                    ind=[regexp(file{1,i},'.ht3') regexp(file{1,i},'.ptu')];
                    if ~isempty(ind)
                        filename{1,i} = [path file{1,i}];
                    elseif LTfile && ~useRaw
                        filename{1,i} = [path file{1,i}];
                    end
                end
            else
                if file~=0
                % raw data is in .ht3-file, user-defined data may have any name
                    ind=[regexp(file,'.ht3') regexp(file,'.ptu')];
                    if ~isempty(ind)
                        filename = [path file];
                    elseif LTfile && ~useRaw
                        filename = [path file]; 
                    end
                end
            end
            if exist('filename','var') % did the user choose a valid file format?
                if objectHandle==app.ChooseFileButton % which button did the user click?
                    app.Filename = filename;
                    set(app.Label_Filename,'String',filename);
                elseif objectHandle==app.IRFFileButton
                    app.IRFFilename = filename;
                    set(app.Label_IRFFilename,'String',filename);
                elseif objectHandle==app.FreeLTButton
                    app.FreeLTFilename = filename;
                    set(app.Label_FreeSpaceLTFilename,'String',filename);
                end
            else
                errordlg('You have to choose an ht3- or a ptu-file.');
            end
        end
        
 
        % %  functions for actual calculations  % % % % % % % % % % % % % %
        function btnCallback(app,~,~) % runs when button is clicked in MIET-modus
        % does app.Filename contain a valid path? does an IRF file exist
        % if we need one?
            if (get(app.MIETEvaluationType,'Value')~=3 && isempty(app.Filename)) || (get(app.IRFcalc,'Value')==1 && isempty(app.IRFFilename))
                errordlg('You have not chosen a valid path to your data.')
                return 
            end
        % cutoff-time
            cutoff = str2double(get(app.CutoffTime,'String'));
        % should the intensity image be made from all photons or only from the photons after the cutoff-time?
            if get(app.CutoffImage,'Value')==1
                intensityImageFlag=1;
            else
                intensityImageFlag=0;
            end
        % if user only wants intensity and lifetime image: do that and be finished
            if get(app.MIETEvaluationType,'Value')==2
                set(app.Figure,'HandleVisibility','off'); % to prevent GUI from being closed by "close all" in subfunctions
            % calculate lifetimes & intensities, are saved in files with
            % the same beginning but ending in '_DATA.mat' and '_tau.mat'
                if ~iscell(app.Filename) % needed to be able to cope both with 1 or several files
                    app.Filename={app.Filename};
                end
                for fileCounter=1:size(app.Filename,2)
                    filename=app.Filename{1,fileCounter};
                    if get(app.LTcalc,'Value')==1 % use raw .ht3-data for calculation
                        [data, tag, life_imm, life_imm_c] = Luminosa_Process_scan(filename,1,[],intensityImageFlag,cutoff,1); %#ok<ASGLU>
                        if size(tag,3)>1 % several detectors?
                            [avg_im, usedPhotons] = app.CombineDetectorLifetimeAverage(tag, life_imm);
                            save([filename(1:end-4) '_PS.mat'],'data','tag','life_imm','life_imm_c','avg_im','usedPhotons');
                        else
                            save([filename(1:end-4) '_PS.mat'],'data','tag','life_imm','life_imm_c');
                        end
                    else
                        load(filename);
                        if ~exist('tag','var') || ~exist('life_imm','var')
                            errordlg(['The file has to contain two variables, "tag" (intensity image)' ...
                                'and "life_imm" (lifetime image), both with same size.']);
                            return;
                        end
                        if size(tag,3)>1 % several detectors?
                            if ~exist('avg_im','var') || ~exist('usedPhotons','var')
                                [avg_im, usedPhotons] = app.CombineDetectorLifetimeAverage(tag, life_imm);
                            end
                        end
                    end
                    app.Lifetime=life_imm;
                    app.Intensity=tag;
                    if size(tag,3)>1 % several detectors?
                        for numChan=1:size(tag,3)
                            figure; imagesc(squeeze(tag(:,:,numChan))); title(sprintf('intensity in detector %i [counts]',numChan)); colorbar;
                            figure; imagesc(squeeze(life_imm(:,:,numChan))); title(sprintf('lifetime in detector %i [ns]',numChan)); colorbar;
                        end
                        figure; imagesc(avg_im); title('lifetime averaged over all detectors'); colorbar;
                        figure; mim(usedPhotons); title(gca,'total number of photons used for calculation'); colorbar;
                        figure; mim(avg_im,usedPhotons); title(gca,'lifetimes weighted with intensity');
                    else % only one detector
                        figure; imagesc(tag); title('intensity image'); colorbar;
                        figure; imagesc(life_imm); title('lifetime image'); colorbar;
                    end
                    set(app.Figure,'HandleVisibility','on');
                end
                return
            end
        % calculate MIET calibration curve
            disp('Calculate MIET calibration curve...');
            app.CalculateCalibrationCurve; 
            calibrationCurve = app.Curve.calibrationCurve; % 1st column=height[nm], 2nd column=lifetime[ns]
            limit_LT = app.Curve.limit_LT; % highest lifetime value in the curve
        % if user only wants MIET calibration curve: you are done!
            if get(app.MIETEvaluationType,'Value')==3
                curveOrientation = 'fast_rotating';
                if isfield(app.Curve,'orientationMode')
                    curveOrientation = app.Curve.orientationMode;
                end
                figure; plot(calibrationCurve(:,2),calibrationCurve(:,1),'.r');
                ylabel('height [nm]'); xlabel('lifetime [ns]');
                title(sprintf('MIET calibration curve (%s)', strrep(curveOrientation, '_', ' ')));
            % re-enable usage of GUI controls (buttons etc.)
                set(app.Figure,'HandleVisibility','on');
                set(app.Button,'enable','on');
                set(app.LTcalc,'enable','on');
                set(app.LTuseOther,'enable','on');
                set(app.ChooseFileButton,'enable','on');
                set(app.MakeImageButton,'enable','on');
                set(app.MakeImages,'enable','on');    
                return;
            end
        % normal calculation
            set(app.Figure,'HandleVisibility','off'); % to prevent GUI from being closed by "close all" in subfunctions
            if ~iscell(app.Filename) % needed to be able to cope both with 1 or several files
                app.Filename={app.Filename};
            end
            for fileCounter=1:size(app.Filename,2)
                filename=app.Filename{1,fileCounter};
                if get(app.LTcalc,'Value')==1 % use raw .ht3-data for calculation
                    if exist([filename(1:end-4) '_PS.mat'],'file')
                    % load data
                        clear data tag life_imm life_imm_c avg_im usedPhotons
                        load([filename(1:end-4) '_PS.mat']); % tag contains intensities, life_imm the lifetimes
                        hasCacheVars = exist('data','var') && exist('tag','var') && exist('life_imm','var');
                        if ~hasCacheVars || app.ProcessScanCacheNeedsRefresh(data, tag, life_imm, cutoff, intensityImageFlag)
                            [data, tag, life_imm, life_imm_c] = Luminosa_Process_scan(filename,1,[],intensityImageFlag,cutoff,1);
                            save([filename(1:end-4) '_PS.mat'],'data','tag','life_imm','life_imm_c');
                        end
                    else
                    % calculate lifetimes & intensities, are saved in file with the
                    % same beginning but ending in '_PS.mat'
                        [data, tag, life_imm, life_imm_c] = Luminosa_Process_scan(filename,1,[],intensityImageFlag,cutoff,1);
                        save([filename(1:end-4) '_PS.mat'],'data','tag','life_imm','life_imm_c');
                    end
                else % use previously calculated values
                    load(filename);
                    if ~exist('tag','var') || ~exist('life_imm','var')
                        errordlg(['The file has to contain two variables, "tag" (intensity image)' ...
                            'and "life_imm" (lifetime image), both with same size.']);
                        return;
                    end
                    if max(life_imm(~isnan(life_imm)))<1e-3 % results were saved in seconds, not nanoseconds
                        life_imm = life_imm*1e9;    % convert to nanoseconds
                    end
                end
                height=zeros(size(life_imm));    % will save result, i.e. height of points
                for i=1:numel(life_imm)
                % Test if lifetime is above unique range or smaller
                % than smallest value of calibration curve
                    if life_imm(i)>limit_LT || life_imm(i)<calibrationCurve(1,2)
                        height(i)=NaN;
                    else
                    % Find 1st element of curve that is larger than current matrix
                    % element, interpolate between it and the previous curve element.
                        for y=2:length(calibrationCurve) 
                            if calibrationCurve(y,2) >= life_imm(i)
                                height(i)=calibrationCurve(y-1,1)+(calibrationCurve(y,1)-calibrationCurve(y-1,1))...
                                    /(calibrationCurve(y,2)-calibrationCurve(y-1,2))*(life_imm(i)-calibrationCurve(y-1,2)); 
                                break; 
                            end
                        end
                    end
                end
                save([filename(1:end-4) '_height.mat'],'height');
                app.Lifetime=life_imm;
                app.Height=height;
                app.Intensity=tag;
            % always make basic plots
                if size(tag,3)>1 % several detectors?
                    avg_im=zeros(size(tag(:,:,1))); avg_height=avg_im; usedPhotons=zeros(size(tag(:,:,1)));
                    for numChan=1:size(tag,3)
                        figure; imagesc(squeeze(tag(:,:,numChan))); title(sprintf('intensity in detector %i [counts]',numChan)); colorbar;
                        figure; imagesc(squeeze(life_imm(:,:,numChan))); title(sprintf('lifetime in detector %i [ns]',numChan)); colorbar;
                        figure; imagesc(squeeze(height(:,:,numChan))); title(sprintf('height of molecules from bottom of their layer in channel %i [nm]',numChan)); colorbar;
                        enoughPhotons=tag(:,:,numChan)>500; eP3d = false(size(tag)); eP3d(:,:,numChan)=enoughPhotons;
                        avg_im(enoughPhotons)  =   avg_im(enoughPhotons)     +tag(eP3d).*life_imm(eP3d);
                        avg_height(enoughPhotons)= avg_height(enoughPhotons) +tag(eP3d).*height(eP3d);
                        usedPhotons(enoughPhotons)=usedPhotons(enoughPhotons)+tag(eP3d);
                    end
                    avg_im = avg_im./usedPhotons;           % image of average lifetimes
                    avg_height = avg_height./usedPhotons;   % image of average heights
                    figure; imagesc(avg_im); title('lifetime averaged over all detectors'); colorbar;
                    figure; imagesc(avg_height); title('height averaged over all detectors'); colorbar;
                    figure; mim(usedPhotons); title(gca,'total number of photons used for calculation'); colorbar;
                    figure; mim(avg_im,usedPhotons); title(gca,'lifetimes weighted with intensity');
                    save([filename(1:end-4) '_averages.mat'],'avg_im','avg_height','usedPhotons');
                else
                   figure; imagesc(tag); title('intensity image'); colorbar;
                   figure; imagesc(life_imm); title('lifetime image'); colorbar;
                   figure; imagesc(height); title('height image, adjust intensity threshold for ignoring dim pixels'); colorbar;
                   uicontrol('Style','edit','String','1','Callback',@app.ThresholdField); % text field for threshold
                end
            % make nice plots if the user wants that
                if ~isempty('app.Height') % did everything work? i.e.: does height exist as a variable?
                    if get(app.MakeImages,'Value')==1 % does the user want to visualize the data?
                    % make sure tag and life_imm have same dimension (binning!)
                        if size(tag,1) ~= size(height,1)
                            binning=floor(size(tag,1)/size(height,1));
                            tag = shiftdim(sum(reshape(tag(1:size(height,1)*binning,:),binning,size(height,1),size(tag,2)),1));
                        end
                        if size(tag,2) ~= size(height,2)
                            binning=floor(size(tag,1)/size(height,1));
                            tag = permute(tag,[2 1 3]);
                            tag = shiftdim(sum(reshape(tag(1:size(height,2)*binning,:),binning,size(height,2),size(tag,2)),1));
                            tag = permute(tag,[2 1 3]);
                        end
                    % find the indices of background pixels, i.e. pixels with a
                    % low intensity: their lifetime is wrong, don't plot
                        thres_index= tag<str2double(get(app.Threshold,'String'));
                        heighttrunc=app.Height;     % height data
                        heighttrunc(thres_index)=NaN; % only plot where thres_index==false
                    % make an array with "shadow" values for a nice 3d plot
                        shadow=ones(size(app.Height))*max(max(heighttrunc))*27/64; % 27/64 will give nice colour 
                        shadow(thres_index)=0;
                    % plot the data
                        makenima(shadow,heighttrunc);
                    end
                end
            end
        % re-enable usage of GUI controls (buttons etc.)
            set(app.Figure,'HandleVisibility','on');
            set(app.Button,'enable','on');
            set(app.LTcalc,'enable','on');
            set(app.LTuseOther,'enable','on');
            set(app.ChooseFileButton,'enable','on');
            set(app.MakeImageButton,'enable','on');
            set(app.MakeImages,'enable','on');
        end
        
        function smMIETCallback(app,~,~) % runs when button is clicked in smMIET-modus
            EvaluationType = get(app.smMIETEvaluationType,'Value');
        % check validity of files (if you need them)           
            if EvaluationType~=8 && ( isempty(app.Filename) || (get(app.IRFcalc,'Value')==1 && isempty(app.IRFFilename)) )
                errordlg('You have not chosen a valid path to your data.')
                return
            elseif EvaluationType~=8
                if ~iscell(app.Filename) % needed to be able to cope both with 1 or several files
                    app.Filename={app.Filename};
                end
            end
        % check if parameters were set (if you need them)
            if EvaluationType~=5 && EvaluationType~=7 && ~isfield(app.Parameters,'NA')
                errordlg('Sample parameters have not been set.');
                return
            elseif EvaluationType~=5 && EvaluationType~=7
            % set parameters; excitation: lengths in µm, emission: lengths in nm
                dye.lamex   = app.Parameters.ExcWavelength*1e-3;
                if app.Parameters.Monochrome
                    dye.lamem   = app.Parameters.EmiWavelength;
                else
                    dye.lamem.Wavel_Small = app.Parameters.Wavel_Small;
                    dye.lamem.Wavel_Large = app.Parameters.Wavel_Large;
                    dye.lamem.SpectrumFile = app.Parameters.SpectrumFile;
                end
                dye.qy      = app.Parameters.QuantumYield;
                dye.tau_free= str2double(get(app.FreeLifetime,'String'));
                dye.CurveType=app.Parameters.CurveType;
                if isfield(app.Parameters,'DipoleOrientation')
                    dye.DipoleOrientation = app.Parameters.DipoleOrientation;
                else
                    dye.DipoleOrientation = 'fast_rotating';
                end

                layers.n1   = app.Parameters.RefrIndexUpper;
                layers.n0   = app.Parameters.RefrIndexLower;
                layers.n    = app.Parameters.RefrIndexMolec;
                layers.d1   = app.Parameters.ThicknessUpper*1e-3;
                layers.d0   = app.Parameters.ThicknessLower*1e-3;
                layers.d    = app.Parameters.ThicknessMolec*1e-3;

                mic.NA      = app.Parameters.NA;
                mic.focpos  = app.Parameters.Focus*1e-3;
                switch        app.Parameters.LaserMode
                    case 1
                        mic.pattern = 'radial';
                    case 2
                        mic.pattern = 'azimuthal';
                    case 3
                        mic.pattern = 'linear';
                end
            end
        % which evaluation type does the user want?
            switch get(app.smMIETEvaluationType,'Value')
                case 1
                    flag = 'MIET';
                case 2
                    flag = 'ROI MIET';
                case 3
                    flag = 'lifetime';
                case 4
                    flag = 'random';
                case 5
                    flag = 'ROI FLIM'; 
                    dye=[]; layers=[]; mic=[]; % no parameters needed
                case 6
                    flag = 'pattern match';
                case 7
                    flag = 'intensity'; % no IRF and no parameters needed
                    for fileCounter=1:size(app.Filename,2)
                        Luminosa_MIET_Analysis(app.Filename{1,fileCounter},[],[],[],[],flag,[],1);
                    end
                    return
                case 8
                    flag = 'show MIET curves'; % no data file and no IRF needed
            end
        % actual calculations
            set(app.Figure,'HandleVisibility','off'); % to prevent GUI from being closed by "close all" in subfunctions
            if get(app.smMIETEvaluationType,'Value')==8 % only show calibration curve
                Luminosa_MIET_Analysis([],[],dye,layers,mic,flag,1); 
            elseif (get(app.IRFcalc,'Value')==1)
                % calculate IRF from extra data
                for fileCounter=1:size(app.Filename,2)
                    Luminosa_MIET_Analysis(app.Filename{1,fileCounter},app.IRFFilename,dye,layers,mic,flag,[],1);
                end
            else
                % estimate IRF from data
                for fileCounter=1:size(app.Filename,2)
                    Luminosa_MIET_Analysis(app.Filename{1,fileCounter},[],dye,layers,mic,flag,[],1);
                end
            end
            set(app.Figure,'HandleVisibility','on');
        end
  
        function CalculateCalibrationCurve(app,~,~) % calculate MIET-curve
        % get setup parameters from the structure app.Parameters
            if ~isfield(app.Parameters,'RefrIndexUpper')
                errordlg('Sample parameters have not been set.');
                return
            end
            n1=app.Parameters.RefrIndexUpper;   % refr. indices above molecule
            d1=app.Parameters.ThicknessUpper;   % thicknesses above molecule
            n0=app.Parameters.RefrIndexLower;   % refr. indices below molecule
            d0=app.Parameters.ThicknessLower;   % thicknesses below molecule
            n =app.Parameters.RefrIndexMolec;   % refr. index of molecule's layer
            d =app.Parameters.ThicknessMolec;   % thickness of molecule's layer
            QY=app.Parameters.QuantumYield;     % free space quantum yield
            orientationMode = 'fast_rotating';
            if isfield(app.Parameters, 'DipoleOrientation')
                orientationMode = luminosa_miet_normalize_orientation(app.Parameters.DipoleOrientation);
            end
            lifetime=str2double(get(app.FreeLifetime,'String'));    % free space lifetime
            if app.Parameters.Monochrome    % monochromatic evaluation
                k=2*pi/app.Parameters.Wavelength;       % vacuum wavevector
                z_start=app.Parameters.Z_Start;         % position of molecule above bottom stack:
                z_stop =app.Parameters.Z_Stop;          % z_numSteps values linearly spaced 
                z_numSteps=app.Parameters.Z_NumSteps;   % from z_start to z_stop
                z=z_start+(z_stop-z_start)*(0:z_numSteps-1)./(z_numSteps-1); % (row vector for LifetimeL)
            else                            % polychromatic evaluation
                wavel_small=app.Parameters.Wavel_Small;
                wavel_large=app.Parameters.Wavel_Large;
                if wavel_large>800 % 'metals.mat' only contains values for lambda=200-800nm
                    disp('Warning: Only wavelengths <800nm are taken into account.');
                    wavel_large=800;
                end
                if wavel_small<200
                    disp('Warning: Only wavelengths >200nm are taken into account.');
                    wavel_small=200;
                end
                if isempty(app.Parameters.SpectrumFile); 
                    errordlg('Please specify the file containing the emission spectrum.');
                    return
                end
                emissionSpectrum=dlmread(app.Parameters.SpectrumFile); % read in raw emission spectrum, ... 
                emissionSpectrum=emissionSpectrum(emissionSpectrum(:,1)>=wavel_small,1:2); % crop to wavelengths
                emissionSpectrum=emissionSpectrum(emissionSpectrum(:,1)<=wavel_large,:);   % fitting the filters
                if isempty(emissionSpectrum)
                    errordlg('The chosen filter settings do not match the emission spectrum.');
                    return
                end
                emissionSpectrum(:,1)=round(emissionSpectrum(:,1)/5);   % if you multiply this by 5, you get wavelengths grouped in 5nm steps
                spectrIntensity=accumarray(emissionSpectrum(:,1)-emissionSpectrum(1,1)+1, emissionSpectrum(:,2)); % sum intensities of grouped wavelengths
                emissionSpectrum=[unique(emissionSpectrum(:,1)*5) spectrIntensity/sum(spectrIntensity)]; % normalise spectrum to get probabilities
            end
        % check if parameters are correct, if yes: start calculation
            if numel(n1)~=numel(d1)+1 || numel(n0)~=numel(d0)+1
                errordlg(['For each stack of materials (below or above the '...
                    'molecule), there has to be one more index of refraction '...
                    'than thickness values. Zeros in the last line are ignored.'])
                return
            end
        % calculate calibration curve: 
        % v = vertical dipole, p = parallel dipole, 
        % d = emission into lower halfspace, u = emission into upper halfspace,
        % all 4 quantities are row vectors
            if app.Parameters.Monochrome
                if strcmp(app.Parameters.CurveType,'maximum') % calculate up to 1st maximum of curve
                    z=0.1:1000;
                    [~,~,~,~,qvd,qvu,qpd,qpu]=LifetimeLSimpsExp(z*k,n0,n,n1,d0*k,d*k,d1*k);
                    calibrationCurve = zeros(numel(z),2);
                    calibrationCurve(:,1) = z;
                    calibrationCurve(:,2) = orientLifetimeCurve(qvu + qvd, qpu + qpd);
                    peak = find(diff(calibrationCurve(:,2))<0,1); % first local maximum of the curve
                    if ~isempty(peak)
                        calibrationCurve = calibrationCurve(1:peak,:);
                    end
                    limit_LT=calibrationCurve(end,2); % largest lifetime
                elseif strcmp(app.Parameters.CurveType,'minimum'); % only use unambiguous lifetime values
                    z=0.1:1000;
                    [~,~,~,~,qvd,qvu,qpd,qpu]=LifetimeLSimpsExp(z*k,n0,n,n1,d0*k,d*k,d1*k);
                    calibrationCurve = zeros(numel(z),2);
                    calibrationCurve(:,1) = z;
                    calibrationCurve(:,2) = orientLifetimeCurve(qvu + qvd, qpu + qpd);
                    peak = find(diff(calibrationCurve(:,2))<0,1); % first local maximum of the curve
                    if isempty(peak)
                        limit_LT = calibrationCurve(end,2);
                    else
                        limit_LT = min(calibrationCurve(peak:end,2)); % lifetime up to which values are unique
                        calibrationCurve = calibrationCurve(1:find(calibrationCurve(:,2)>limit_LT,1)-1,:);
                    end
                else % use height values chosen by user
                    [~,~,~,~,qvd,qvu,qpd,qpu]=LifetimeLSimpsExp(z*k,n0,n,n1,d0*k,d*k,d1*k);
                    calibrationCurve = zeros(numel(z),2);
                    calibrationCurve(:,1) = z;
                    calibrationCurve(:,2) = orientLifetimeCurve(qvu + qvd, qpu + qpd);
                    % lifetime grows monotonously, then oscillates: 
                    % find height value up to which the mapping lifetime->height is unique
                    peak = find(diff(calibrationCurve(:,2))<0,1); % starting point of oscillations
                    limit_LT = calibrationCurve(end,2); % if no oscillations: last LT is largest one
                    if ~isempty(peak) % do we even simulate up to the oscillations?
                        limit_LT = min(calibrationCurve(peak:end,2)); % lifetime up to which values are unique
                    end
                end
            else % polychromatic evaluation
                z=0.1:1000; % distance between molecule and interface
                calibrationCurve = zeros(numel(z),2); % container to store average calibration curve
                calibrationCurve(:,1) = z;
                n0_backup=n0; n1_backup=n1;  % refractive indices are overwritten for different wavelengths
                reverseStr = '';
                for WVL_index=1:size(emissionSpectrum,1) % loop through all wavelengths
                    lambda = emissionSpectrum(WVL_index,1); % vacuum wavelength
                    k=2*pi/lambda;                          % vacuum wavevector
                    for n_index=1:numel(n0) % replace placeholders by correct refr. indices
                       if(n0_backup(n_index)>=10)
                          n0(n_index)=refrIndex(n0_backup(n_index)); 
                       end
                    end
                    for n_index=1:numel(n1) % replace placeholders by correct refr. indices
                       if(n1_backup(n_index)>=10)
                          n1(n_index)=refrIndex(n1_backup(n_index)); 
                       end
                    end
                    [~,~,~,~,qvd,qvu,qpd,qpu]=LifetimeLSimpsExp(z*k,n0,n,n1,d0*k,d*k,d1*k);
                    calibrationCurve(:,2) = calibrationCurve(:,2) + emissionSpectrum(WVL_index,2)...
                        * orientLifetimeCurve(qvu + qvd, qpu + qpd); % weighted average of calibr. curves
                    msg = sprintf('Processed %d/%d wavelengths...\n', WVL_index, size(emissionSpectrum,1));
                    fprintf([reverseStr, msg]);
                    reverseStr = repmat(sprintf('\b'), 1, length(msg));
                end
                if strcmp(app.Parameters.CurveType,'maximum') % calculate up to 1st maximum of curve
                    peak = find(diff(calibrationCurve(:,2))<0,1); % first local maximum of the curve
                    if ~isempty(peak)
                        calibrationCurve = calibrationCurve(1:peak,:);
                    end
                    limit_LT = calibrationCurve(end,2);           % largest lifetime in calibration curve
                elseif strcmp(app.Parameters.CurveType,'minimum'); % only use unambiguous lifetime values
                    peak = find(diff(calibrationCurve(:,2))<0,1); % first local maximum of the curve
                    if isempty(peak)
                        limit_LT = calibrationCurve(end,2);
                    else
                        limit_LT = min(calibrationCurve(peak:end,2)); % lifetime up to which values are unique
                        calibrationCurve = calibrationCurve(1:find(calibrationCurve(:,2)>limit_LT,1)-1,:);
                    end
                else % use user-defined interval (BUT NOT USING NUMSTEPS STEPS, DIFFERENCE IS ALWAYS 1nm!)
                    calibrationCurve = calibrationCurve(calibrationCurve(:,1)>=app.Parameters.Z_Start,:);
                    calibrationCurve = calibrationCurve(calibrationCurve(:,1)<=app.Parameters.Z_Stop,:);
                    peak = find(diff(calibrationCurve(:,2))<0,1); % starting point of oscillations
                    limit_LT = calibrationCurve(end,2); % if no oscillations: last LT is largest one
                    if ~isempty(peak) % do we even simulate up to the oscillations?
                        limit_LT = min(calibrationCurve(peak:end,2)); % lifetime up to which values are unique
                    end
                end
            end
        % save and display the calibration curve
            save('calibrationCurve.mat','calibrationCurve');
            disp('Calibration curve:');
            disp(calibrationCurve);
            app.Curve.calibrationCurve = calibrationCurve;
            app.Curve.limit_LT = limit_LT;
            app.Curve.orientationMode = orientationMode;
            
            function value = refrIndex(material) % nested function: find refractive index of "material" at wavelength lambda
                switch material
                    case 10
                        value=app.Metals.silver(app.Metals.wavelength==lambda);
                    case 20
                        value=app.Metals.gold(app.Metals.wavelength==lambda);
                    case 30
                        value=app.Metals.platinum(app.Metals.wavelength==lambda);
                    case 40
                        value=app.Metals.palladium(app.Metals.wavelength==lambda);
                    case 50
                        value=app.Metals.copper(app.Metals.wavelength==lambda);
                    case 60
                        value=app.Metals.aluminum(app.Metals.wavelength==lambda);
                    case 70
                        value=app.Metals.chromium(app.Metals.wavelength==lambda);
                    case 80
                        value=app.Metals.titan(app.Metals.wavelength==lambda);
                    case 90
                        value=app.Metals.tungsten(app.Metals.wavelength==lambda);
                    case 100
                        value=app.Metals.nickel(app.Metals.wavelength==lambda);
                    case 110
                        value=app.Metals.beryllium(app.Metals.wavelength==lambda);
                    case 120
                        value=app.Metals.ito(app.Metals.wavelength==lambda);
                    otherwise
                        value=[]; 
                        errormsg('The chosen material does not exist.');
                end
            end

            function values = orientLifetimeCurve(sv, sp)
                switch orientationMode
                    case 'vertical'
                        sr = sv(:);
                        values = lifetime ./ ((1 - QY) + QY * sr / (4/3 * n));
                    case 'parallel'
                        sr = sp(:);
                        values = lifetime ./ ((1 - QY) + QY * sr / (4/3 * n));
                    case 'random_fixed'
                        thetaGrid = linspace(0, pi/2, 181);
                        weights = sin(thetaGrid);
                        sr = sv(:) * (cos(thetaGrid).^2) + sp(:) * (sin(thetaGrid).^2);
                        lifeByTheta = lifetime ./ ((1 - QY) + QY * sr / (4/3 * n));
                        values = trapz(thetaGrid, lifeByTheta .* weights, 2) ./ trapz(thetaGrid, weights);
                    otherwise
                        values = lifetime ./ (1 - QY + QY * (sv(:) + 2 .* sp(:)) / (4 * n));
                end
            end
        end
        
        
        function MakeImageBtnCallback(app,~,~) % plot data manually
        % do app.Height and app.Intensity contain values?
            if ~isempty(app.Height) && ~isempty(app.Intensity) 
                tag=app.Intensity; tav=app.Height;
        % use lifetime and intensity instead
            elseif ~isempty(app.Lifetime) && ~isempty(app.Intensity)
                tag=app.Intensity; tav=app.Lifetime;
            else
                warndlg('The converted data or the intensity values are missing.');
            end
        % make sure Height and Intensity have same dimension (binning!)
            if size(tag,1) ~= size(tav,1)
                binning=floor(size(tag,1)/size(tav,1));
                tag = shiftdim(sum(reshape(tag(1:size(tav,1)*binning,:),binning,size(tav,1),size(tag,2)),1));
            end
            if size(tag,2) ~= size(tav,2)
                binning=floor(size(tag,1)/size(tav,1));
                tag = permute(tag,[2 1 3]);
                tag = shiftdim(sum(reshape(tag(1:size(tav,2)*binning,:),binning,size(tav,2),size(tav,1)),1));
                tag = permute(tag,[2 1 3]);
            end
            app.Intensity=tag;
        % find the indices of background pixels, i.e. pixels with a
        % low intensity: their lifetime is wrong, don't plot
            thres_index= app.Intensity<str2double(get(app.Threshold,'String'));
            heighttrunc=tav;     % height or lifetime data
            heighttrunc(thres_index)=NaN; % only plot where thres_index==false
        % make an array with "shadow" values for a nice 3d plot
            shadow=ones(size(tav))*max(max(heighttrunc))*27/64; % 27/64 will give nice colour 
            shadow(thres_index)=0;
        % plot the data
            makenima(shadow,heighttrunc);
        end
        
        function GetFreeSpaceLifetime(app,~,~) % display estimates of free space lifetime
        % does app.FreeLTFilename contain a valid path?
            if isempty(app.FreeLTFilename)
                errordlg('You have not chosen a valid path to your data.')
            end
        % get free space lifetime from tcspc histogram of the whole image (tau_fit, tau_avg)
        % and from tcspc histograms of single pixels (tau_pixels)
            [tau_fit,tau_avg,tau_px_max,tau_px_mean] = FreeSpaceLifetime(app.FreeLTFilename); 
            fprintf('\ntau_fit = %.2f; tau_avg = %.2f; tau_px_max = %.2f; tau_px_mean = %.2f\n\n',...
                tau_fit,tau_avg,tau_px_max,tau_px_mean);
        end
        
        function ThresholdField(app,ObjectHandle,~) % do not show dim pixels in height image
            threshold = str2double(get(ObjectHandle,'String'));
            tmp = app.Height; tmp(app.Intensity<threshold)=NaN;
            imagesc(tmp); colorbar;
            title('height image, adjust intensity threshold for ignoring dim pixels');
        end

        function [avg_im, usedPhotons] = CombineDetectorLifetimeAverage(app, tag, life_imm) %#ok<INUSD>
            avg_im = zeros(size(tag(:,:,1)));
            usedPhotons = zeros(size(tag(:,:,1)));
            for numChan = 1:size(tag,3)
                enoughPhotons = tag(:,:,numChan) > 500;
                eP3d = false(size(tag));
                eP3d(:,:,numChan) = enoughPhotons;
                avg_im(enoughPhotons) = avg_im(enoughPhotons) + tag(eP3d) .* life_imm(eP3d);
                usedPhotons(enoughPhotons) = usedPhotons(enoughPhotons) + tag(eP3d);
            end
            valid = usedPhotons > 0;
            avg_im(valid) = avg_im(valid) ./ usedPhotons(valid);
            avg_im(~valid) = NaN;
        end

        function tf = ProcessScanCacheNeedsRefresh(app, data, tag, life_imm, cutoff, intensityImageFlag) %#ok<INUSD>
            tf = ~isstruct(data) || ~isfield(data,'cutoff') || data.cutoff ~= cutoff || ...
                ~isfield(data,'intensityImageFlag') || data.intensityImageFlag ~= intensityImageFlag || ...
                ~isfield(data,'lifetimeEstimator') || ~strcmp(data.lifetimeEstimator, app.LocalProcessScanEstimatorLabel()) || ...
                isempty(tag) || isempty(life_imm);
        end

        function estimatorLabel = LocalProcessScanEstimatorLabel(app) %#ok<MANU>
            estimatorLabel = 'Luminosa_DistFluofit_extension_PIRLSGPU_v1';
        end
    end
      
end
