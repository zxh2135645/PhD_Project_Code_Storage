classdef GUI2022 < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                    matlab.ui.Figure
        GridLayout                  matlab.ui.container.GridLayout
        LeftPanel                   matlab.ui.container.Panel
        ImageModalityDropDownLabel  matlab.ui.control.Label
        ImageModalityDropDown       matlab.ui.control.DropDown
        NameListDropDownLabel       matlab.ui.control.Label
        NameListDropDown            matlab.ui.control.DropDown
        LoadFileButton              matlab.ui.control.Button
        SaveFileButton              matlab.ui.control.Button
        LoadCorrectedMaskButton     matlab.ui.control.Button
        SliceEditFieldLabel         matlab.ui.control.Label
        SliceEditField              matlab.ui.control.EditField
        EnhancedorNotButtonGroup    matlab.ui.container.ButtonGroup
        NotButton                   matlab.ui.control.RadioButton
        EnhancedButton              matlab.ui.control.RadioButton
        SliceskippedEditFieldLabel  matlab.ui.control.Label
        SliceskippedEditField       matlab.ui.control.EditField
        ImageSavedStatusLampLabel   matlab.ui.control.Label
        ImageSavedStatusLamp        matlab.ui.control.Lamp
        WindowingPanel              matlab.ui.container.Panel
        UpperSliderLabel            matlab.ui.control.Label
        UpperSlider                 matlab.ui.control.Slider
        LowerSliderLabel            matlab.ui.control.Label
        LowerSlider                 matlab.ui.control.Slider
        RightPanel                  matlab.ui.container.Panel
        OrigImagePanel              matlab.ui.container.Panel
        EditPanel                   matlab.ui.container.Panel
        %RectRemoveButton            matlab.ui.control.Button
        PixelRemoveButton           matlab.ui.control.StateButton
        ROIRemoveButton             matlab.ui.control.Button
        %RectAddButton               matlab.ui.control.Button
        PixelAddButton              matlab.ui.control.StateButton
        ROIAddButton                matlab.ui.control.Button
        RedoButton                  matlab.ui.control.Button
        SkipButton                  matlab.ui.control.Button
        ModeButtonGroup             matlab.ui.container.ButtonGroup
        EndoButton                  matlab.ui.control.ToggleButton
        EpiButton                   matlab.ui.control.ToggleButton
        ExclusionButton             matlab.ui.control.ToggleButton
        RemoteButton                matlab.ui.control.ToggleButton
        ReflowButton                matlab.ui.control.ToggleButton
        Panel                       matlab.ui.container.Panel
        Panel_2                     matlab.ui.container.Panel
        PreviousButton_2            matlab.ui.control.Button
        NextButton_2                matlab.ui.control.Button
        AlgorithmsDropDownLabel     matlab.ui.control.Label
        AlgorithmsDropDown          matlab.ui.control.DropDown
        ComputeButton               matlab.ui.control.Button
        Slider                      matlab.ui.control.Slider
        ImageOverviewPanel          matlab.ui.container.Panel
        TabGroup                    matlab.ui.container.TabGroup
        Tab                         matlab.ui.container.Tab
        Tab2                        matlab.ui.container.Tab
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
        funcPath
        % Input parameters
        % They are hard coded, so if this app is for other purposes, these
        % are needed to be taken into consideration carefully.
        % sequence_label = {'D1', 'D2', 'D3'};
        image_modal
        popu = 'Patient';
        enh = 'Not';
        name_glob
        tp_glob
        name
        tp
        real_pwd
        base_dir
        label
        dicom_glob
        slice_info
        sloc_array
        subject_idx
        
        % Initialize before load
        dicom_fields = {...
            'Filename',...
            'Height', ...
            'Width', ...
            'Rows',...
            'Columns', ...
            'PixelSpacing',...
            'SliceThickness',...
            'SliceLocation',...
            'ImagePositionPatient',...
            'ImageOrientationPatient',...
            'MediaStorageSOPInstanceUID',...
            };
        
        img_3D
        img_4D
        img_size
        slice_idx
        rgb_scale
        ax1_rgb
        ax1
        ax1_2
        img_3D_rgb
        ax2
        ax3
        
        mode_btn = 'Endo';
        
        % Edit functions
        window_button_down
        rect_mask_button_down
        rect_mask_add_button_down
        rect_or_pixel
        loc_idx
        loc
            
        % Rect and pixel
        rect_mask_count
        rect_mask_add_count
        pixel_mask_rm_count
        pixel_mask_add_count
        
        % ROI functions
        roi_mask_add_count
        roi_mask_add
        roi_mask_rm_count
        roi_mask_rm
        
        % Endo, Epi, Exclusion, Remote
        endo_mask
        endo_mask_moded
        epi_mask
        epi_mask_moded
        exclusion_mask
        exclusion_mask_moded
        remote_mask
        remote_mask_moded
        reflow_mask
        reflow_mask_moded
        
        myo_mask
        myo_mask_moded
        mi_mask
        mi_mask_moded
        composite_mask
        
        first_edit_end_label
        first_edit_epi_label
        first_edit_exc_label 
        first_edit_rem_label
        first_edit_ref_label
        % check if it's first time clicked edit panel per each slice
        % 0 - uneditted    1 - first time edit    2... - 2nd or more edit
        
        % Algorithms label
        algo_label = {'Mean+5SD', 'Mean-2SD', 'Manual'};
        algo
        manual_thresh
        computebtn
        
        % SLICE Skip
        skip_slice_idx
        
        % Output variables
        save_dir
        name_dir
        blood_dir
        exclude_dir
        heart_dir
        myo_dir
        remote_dir
        mi_dir
        myo_pool_dir
        reflow_dir
        
        % Windowing slider
        upper_window
        lower_window
        img_for_window
        
        % Image overview
        ax_stream = struct;
    end

    % Callbacks that handle component events
    methods (Access = private)
        
        function startupFcn(app)
            %if (~isdeployed)
            %    app.funcPath = addpath('functions/'); % Addpath to functions
            %end
            app.Slider.Enable = 'off';
        end
        
        % Value changed function: ImageModalityDropDown
        function ImageModalityDropDownValueChanged(app, event)
            app.name = app.ImageModalityDropDown.Value;
        end

        % Selection changed function: EnhancedorNotButtonGroup
        function EnhancedorNotButtonGroupSelectionChanged(app, event)
            app.enh = app.EnhancedorNotButtonGroup.SelectedObject;
        end

        % Value changed function: NameListDropDown
        function NameListDropDownValueChanged(app, event)
            app.tp = app.NameListDropDown.Value;
        end
        
        function GetSliceInfo(app, event)
            app.slice_info = zeros(1, length(app.dicom_glob));
            for i = 1:length(app.dicom_glob)
                strings = strsplit(app.dicom_glob{i}, '/');
                f = strings{end};
                strings2 = strsplit(f, '_');
                idx_str = strings2{end};
                idx = sscanf(idx_str, '%d*');
                app.slice_info(i) = idx;
            end
        end
        
        function ReadDICOM_3D(app, event)
            app.slice_info = zeros(1, length(app.dicom_glob));
            
            for i = 1:length(app.dicom_glob)
                [X, map] = dicomread(app.dicom_glob{i});
                if i == 1
                   app.img_3D = zeros(size(X,1), size(X,2), length(app.dicom_glob)); 
                end
                strings = strsplit(app.dicom_glob{i}, '/');
                f = strings{end};
                strings2 = strsplit(f, '_');
                idx_str = strings2{end};
                idx = sscanf(idx_str, 'slice%d*');
                app.slice_info(i) = idx;
                
                app.img_3D(:,:,idx) = X;
            end
        end
        
        function GetSubjectIndex(app, event)
            strings = strsplit(app.name, '_');
            f = strings{end};
            app.subject_idx = str2num(f);
        end
        
        function GetRealPWD(app, event)
            if isdeployed && ismac
                NameOfDeployedApp = 'GUI2022'; % do not include the '.app' extension
                [~, result] = system(['top -n100 -l1 | grep ' NameOfDeployedApp ' | awk ''{print $1}''']);
                result=strtrim(result);
                [status, result] = system(['ps xuwww -p ' result ' | tail -n1 | awk ''{print $NF}''']);
                if status==0
                    diridx=strfind(result,[NameOfDeployedApp '.app']);
                    app.real_pwd=result(1:diridx-2);
                else
                    msgbox({'realpwd not set:',result})
                end
            else
                app.real_pwd = pwd;
            end
        end
        
        function T2starOverlayFig(app, event)
            imagesc(app.img_3D_rgb); pbaspect([size(app.img_3D_rgb,2), size(app.img_3D_rgb,1) 1]);
            colormap(app.ax1,'gray')
            set(app.ax1,'xticklabel',[]); set(app.ax1,'yticklabel',[])
            imagesc(app.ax1_2, squeeze(app.myo_mask(:,:,app.slice_idx)), 'AlphaData', squeeze(app.myo_mask(:,:,app.slice_idx))); 
            pbaspect([size(app.img_3D_rgb,2), size(app.img_3D_rgb,1) 1]);
            colormap(app.ax1_2,'jet')
            %caxis(ax1,[0 800])
            %caxis(ax2,[0 100])
            linkprop([app.ax1, app.ax2], 'Position')
            app.ax2.Visible = 'off';
            cb = colorbar;
            title(cb, 'T_{2}* (ms)')
        end
        
        % Button pushed function: LoadFileButton
        function LoadFileButtonPushed(app, event)
            
            GetRealPWD(app);
            app.base_dir = cat(2,  app.real_pwd, '/PatientFatAnalysis/');
            
            % dicom_to_read = cat(2, app.base_dir, 'Dicom/', app.name, '/', app.image_modal, '/');
            f_to_read = cat(2, app.base_dir, 'ContourData/', app.name, '/', app.tp, '/T2star/');
            %app.dicom_glob = glob(cat(2, dicom_to_read, '*.dcm'));
            %GetSliceInfo(app); % Not further used yet
            D = load(cat(2, f_to_read, 'T2star_vol_img_4D.mat'));
            Dx = D.('T2star_struct');
            %[app.img_3D, slice_data, image_meta_data] = dicom23D(app.dicom_glob{1}, app.dicom_fields);
            %ReadDICOM_3D(app);
            
            %app.sloc_array = app.slice_info;
%             for i = 1:length(slice_data)
%                app.sloc_array(i) = slice_data(i).SliceLocation;
%             end
            GetSubjectIndex(app);
            % app.sloc_array = size(Dx{app.subject_idx}, 1);
            app.img_4D = Dx.('vol_img_4D'); % TODO: a 4D t2* weighted image will be implemented
            
            app.img_3D = app.img_4D(:,:,:,1); % It is now hard-coded here
            
            app.slice_idx = 1;
            app.skip_slice_idx = [];
            app.SliceskippedEditField.Value = mat2str(app.skip_slice_idx);
            app.rgb_scale = 0.8;
            app.img_size = [size(app.img_4D,1), size(app.img_4D,2), size(app.img_4D,3)];
            
            % Initialize masks
            app.endo_mask = zeros(app.img_size);
            app.endo_mask_moded = zeros(app.img_size);
            app.epi_mask = zeros(app.img_size);
            app.epi_mask_moded = zeros(app.img_size);
            app.exclusion_mask = zeros(app.img_size);
            app.exclusion_mask_moded = zeros(app.img_size);
            app.remote_mask = zeros(app.img_size);
            app.remote_mask_moded = zeros(app.img_size);
            app.reflow_mask = zeros(app.img_size);
            app.reflow_mask_moded = zeros(app.img_size);
            
            app.myo_mask = app.epi_mask_moded - app.endo_mask_moded;
            app.mi_mask = zeros(app.img_size);
            app.mi_mask_moded = app.mi_mask;
            app.composite_mask = zeros(app.img_size);
            app.myo_mask_moded = app.myo_mask;
            
            app.ax1_rgb = zeros(app.img_size(1), app.img_size(2), 3);
            % app.img_3D_rgb = adapthisteq(app.img_3D(:,:,app.slice_idx) / max(max(app.img_3D(:,:,app.slice_idx))));
            app.img_3D_rgb = app.img_3D(:,:,app.slice_idx) / max(max(app.img_3D(:,:,app.slice_idx)));
            
            app.ax1_rgb(:,:,1) = app.rgb_scale*app.img_3D_rgb + app.mi_mask(:,:,app.slice_idx);
            app.ax1_rgb(:,:,2) = app.rgb_scale*app.img_3D_rgb + app.remote_mask(:,:,app.slice_idx);
            % edg  = edge(squeeze(app.myo_mask(:,:,app.slice_idx)),'Canny');
            app.ax1_rgb(:,:,3) = app.rgb_scale*app.img_3D_rgb + app.myo_mask(:,:,app.slice_idx);
            
            app.ax1 = uiaxes(app.OrigImagePanel,'Position',[58,2,375,390]);
            imagesc(app.ax1_rgb, 'Parent', app.ax1);
            set(app.ax1, 'visible', 'off');
            axis(app.ax1, 'image');
%             
            %imcontrast(app.ax1_rgb);
            
            
            %app.img_3D_rgb = adapthisteq(app.img_3D(:,:,app.slice_idx) / max(max(app.img_3D(:,:,app.slice_idx))));
            %app.ax1 = uiaxes(app.OrigImagePanel,'Position',[58,2,375,390]);
            %app.ax1_2 = uiaxes(app.OrigImagePanel,'Position',[58,2,375,390]);
            %T2starOverlayFig(app);
            
            app.ax2 = uiaxes(app.RightPanel,'Position',[10 1 237 231]);
            imagesc(app.img_3D_rgb, 'Parent', app.ax2);
            set(app.ax2, 'visible', 'off');
            axis(app.ax2, 'image');
            colormap(app.ax2,gray);
            
            app.ax3 = uiaxes(app.RightPanel,'Position',[250 1 237 231]);
            imagesc(app.composite_mask(:,:,app.slice_idx), 'Parent', app.ax3);
            set(app.ax3, 'visible', 'off');
            axis(app.ax3, 'image');
            
            vert_interval = 160;
            for st = 1:app.img_size(3)
                ax_overview = app.img_3D(:,:,st) / max(max(app.img_3D(:,:,st)));
                if st <= 4
                    app.ax_stream(st).ax = uiaxes(app.Tab, 'Position', [10, 469 - vert_interval*(st-1), 212, 200]);
                    imagesc(ax_overview, 'Parent', app.ax_stream(st).ax);
                    set(app.ax_stream(st).ax, 'visible', 'off');
                    colormap(app.ax_stream(st).ax, gray);
                    % axis(app.ax_stream(1), 'image');
                    % title(app.ax_stream(1), 'Slice 1');
                    app.ax_stream(st).label = uilabel(app.Tab);
                    app.ax_stream(st).label.Position = [190, 560-vert_interval*(st-1), 37, 26];
                    app.ax_stream(st).label.Text = num2str(st);
                    app.ax_stream(st).label.FontSize = 18;
                elseif st > 4
                    app.ax_stream(st).ax = uiaxes(app.Tab2, 'Position', [10, 469 - vert_interval*(st-5), 212, 200]);
                    imagesc(ax_overview, 'Parent', app.ax_stream(st).ax);
                    set(app.ax_stream(st).ax, 'visible', 'off');
                    colormap(app.ax_stream(st).ax, gray);
                    app.ax_stream(st).label = uilabel(app.Tab2);
                    app.ax_stream(st).label.Position = [190, 560-vert_interval*(st-5), 37, 26];
                    app.ax_stream(st).label.Text = num2str(st);
                    app.ax_stream(st).label.FontSize = 18;
                end
                
                disableDefaultInteractivity(app.ax_stream(st).ax);
            end
            
            app.Tab.Scrollable = 'on';
            
            % Try to accelarate the program
            disableDefaultInteractivity(app.ax2);
            disableDefaultInteractivity(app.ax3);
            
            % Slice
            slice_idx_str = num2str(app.slice_idx);
            max_slice_idx_str = num2str(app.img_size(3));
            app.SliceEditField.Value = cat(2, slice_idx_str, ' / ', max_slice_idx_str);
            % Image Saved Status
            app.ImageSavedStatusLamp.Color = [1 0 0];
            
            % To make buttons mutually exclusive
            app.window_button_down = 0;
            app.rect_mask_button_down = 0;
            app.rect_mask_add_button_down = 0;
            app.rect_or_pixel = 0;
            
            app.rect_mask_count = 0; % Reset rect_mask_count
            app.rect_mask_add_count = 0;
            app.pixel_mask_rm_count = 0;
            app.pixel_mask_add_count = 0;
            app.roi_mask_rm_count = 0;
            app.roi_mask_add_count = 0;
            
            app.roi_mask_rm = zeros(app.img_size(1), app.img_size(2));
            app.roi_mask_add = zeros(app.img_size(1), app.img_size(2));
            
            % Compute button initialize
            app.computebtn = 0;
            app.manual_thresh = zeros(1, app.img_size(3));
            
            app.UpperSlider.Limits = [min(app.img_3D_rgb(:)) max(app.img_3D_rgb(:))];
            app.LowerSlider.Limits = [min(app.img_3D_rgb(:)) max(app.img_3D_rgb(:))];
            app.UpperSlider.Value = max(app.img_3D_rgb(:));
            app.LowerSlider.Value = min(app.img_3D_rgb(:));
            
            app.upper_window = app.UpperSlider.Value;
            app.lower_window = app.LowerSlider.Value;
            app.img_for_window = app.img_3D_rgb;
        end

        % Button pushed function: SaveFileButton
        function SaveFileButtonPushed(app, event)
            contour_dir = cat(2, app.base_dir, 'ContourData/');
            if ~exist(contour_dir, 'dir')
               mkdir(contour_dir); 
            end
            
            %popu_dir = cat(2, contour_dir, app.popu, '/');
            %if ~exist(popu_dir, 'dir')
            %   mkdir(popu_dir);
            %end
            
            app.name_dir = cat(2, contour_dir, app.name, '/');
            if ~exist(app.name_dir, 'dir')
               mkdir(app.name_dir);
            end
            
            labelo = app.tp;
            app.save_dir = cat(2, app.name_dir, labelo, '/');
            if ~exist(app.save_dir, 'dir')
               mkdir(app.save_dir); 
            end
            
            mi_mask_check = reshape((app.mi_mask_moded), [], app.img_size(3));
            idx = find(any(mi_mask_check) == 0); % idx is any slices that is not edited
            isSubset = all(ismember(idx, app.skip_slice_idx)); % Try to see if not-edited is a subset of skipped
            notMember = setdiff(idx, app.skip_slice_idx);
            
            if all(any(mi_mask_check)) || isSubset % check if all are edited, or if not-edited are skipped
                % Yes mkdir in save_dir
                mkdir_in_save_dir(app);
                mask_blood_3D = app.endo_mask_moded;
                excludeMask_3D = app.exclusion_mask_moded;
                mask_heart_3D = app.epi_mask_moded;
                mask_myocardium_3D = app.myo_mask;
                myoRefMask_3D = app.remote_mask;
                infarct_masked_3D = app.mi_mask_moded;
                reflowMask_3D = app.reflow_mask_moded;
                compositeIm = app.composite_mask;
                vi3 = struct;
                vi3.vol_img_3D = app.img_3D;
                vi3.sliceloc_array = app.sloc_array;
                
                save(cat(2, app.blood_dir, 'mask_blood.mat'), 'mask_blood_3D');
                save(cat(2, app.exclude_dir, 'excludeArea.mat'), 'excludeMask_3D');
                save(cat(2, app.heart_dir, 'mask_heart.mat'), 'mask_heart_3D');
                save(cat(2, app.myo_dir, 'mask_myocardium.mat'), 'mask_myocardium_3D');
                save(cat(2, app.myo_pool_dir, 'mask_myocardium.mat'), 'mask_myocardium_3D');
                save(cat(2, app.remote_dir, 'myoRef.mat'), 'myoRefMask_3D');
                save(cat(2, app.mi_dir, 'Infarct_Whole.mat'), 'infarct_masked_3D');
                save(cat(2, app.reflow_dir, 'noReflowArea.mat'), 'reflowMask_3D');
                save(cat(2, app.save_dir, 'compositeLabel.mat'), 'compositeIm');
                save(cat(2, app.save_dir, app.tp, '_vol_img_3D.mat'), 'vi3');
                % Status Change               
               app.ImageSavedStatusLamp.Color = [0 1 0];
            else
                fig = uifigure;
                message = sprintf('Missing slice at %s.', mat2str(notMember));
                uialert(fig, message, 'Warning', 'Icon', 'warning');
            end
        end

        % Button pushed function: LoadCorrectedMaskButton
        function LoadCorrectedMaskButtonPushed(app, event)
            contour_dir = cat(2, app.base_dir, 'ContourData/');
            app.name_dir = cat(2, contour_dir, app.name, '/');
            load(cat(2, app.name_dir, app.tp, '/T2star/Myocardium/mask_myocardium.mat'));
            app.myo_mask = mask_myocardium_3D;

            % edg  = edge(squeeze(app.myo_mask(:,:,app.slice_idx)),'Canny');
            app.ax1_rgb(:,:,3) = app.rgb_scale*app.img_for_window + app.myo_mask(:,:,app.slice_idx);
            
            update_composite_mask(app);
            refreshplot(app);
            
            check_Manual_Compute(app);
        end

        % Button pushed function: PreviousButton_2
        function PreviousButton_2Pushed(app, event)
             if app.slice_idx > 1
                app.slice_idx = app.slice_idx - 1;
                slice_idx_str = num2str(app.slice_idx);
                max_slice_idx_str = num2str(app.img_size(3));
                app.SliceEditField.Value = cat(2, slice_idx_str, ' / ', max_slice_idx_str);
                
                % RGB update
                % app.img_3D_rgb = adapthisteq(app.img_3D(:,:,app.slice_idx) / max(max(app.img_3D(:,:,app.slice_idx))));
                app.img_3D_rgb = app.img_3D(:,:,app.slice_idx) / max(max(app.img_3D(:,:,app.slice_idx)));
                app.img_for_window = app.img_3D_rgb;
                
                app.UpperSlider.Limits = [min(app.img_3D_rgb(:)) max(app.img_3D_rgb(:))];
                app.LowerSlider.Limits = [min(app.img_3D_rgb(:)) max(app.img_3D_rgb(:))];
                app.UpperSlider.Value = max(app.img_3D_rgb(:));
                app.LowerSlider.Value = min(app.img_3D_rgb(:));
                
                app.upper_window = app.UpperSlider.Value;
                app.lower_window = app.LowerSlider.Value;
                
                app.Slider.Enable = 'off';
                app.Slider.Value = app.manual_thresh(app.slice_idx);
                
                update_all_axes(app);
                
                refreshplot(app);
                
                update_windowing(app);

                check_Manual_Compute(app);
             end
        end

        % Button pushed function: NextButton_2
        function NextButton_2Pushed(app, event)
            if app.slice_idx < app.img_size(3)
                app.slice_idx = app.slice_idx + 1;
                slice_idx_str = num2str(app.slice_idx);
                max_slice_idx_str = num2str(app.img_size(3));
                app.SliceEditField.Value = cat(2, slice_idx_str, ' / ', max_slice_idx_str);
                
                % RGB
                % app.img_3D_rgb = adapthisteq(app.img_3D(:,:,app.slice_idx) / max(max(app.img_3D(:,:,app.slice_idx))));
                app.img_3D_rgb = app.img_3D(:,:,app.slice_idx) / max(max(app.img_3D(:,:,app.slice_idx)));
                app.img_for_window = app.img_3D_rgb;
                
                app.UpperSlider.Limits = [min(app.img_3D_rgb(:)) max(app.img_3D_rgb(:))];
                app.LowerSlider.Limits = [min(app.img_3D_rgb(:)) max(app.img_3D_rgb(:))];
                app.UpperSlider.Value = max(app.img_3D_rgb(:));
                app.LowerSlider.Value = min(app.img_3D_rgb(:));
                
                app.upper_window = app.UpperSlider.Value;
                app.lower_window = app.LowerSlider.Value;
                
                app.Slider.Enable = 'off';
                app.Slider.Value = app.manual_thresh(app.slice_idx);
                
                update_all_axes(app);
                
                refreshplot(app);
                
                update_windowing(app);
                
                check_Manual_Compute(app);
            end
        end
        
        function UIFigureWindowButtonDown(app, event)
            % Need to check if this is called
            if app.rect_mask_button_down == 1
                app.window_button_down = 1;
            elseif app.rect_mask_add_button_down == 1
                app.window_button_down = 2;
            elseif app.PixelRemoveButton.Value == 1
                app.window_button_down = 3;
            elseif app.PixelAddButton.Value == 1
                app.window_button_down = 4;
            end
            
            if app.rect_or_pixel == 1
                % two points for rect
                if app.loc_idx < 2
                    % Fourth quadrant need to pay attention to that
                    app.loc_idx = app.loc_idx + 1;
                    hold(app.ax1, 'on')
                    temp = app.ax1.CurrentPoint; % Returns 2x3 array of points
                    app.loc(:,app.loc_idx) = [temp(1,1); temp(1,2)]; % Gets the (x,y) coordinates
                    plot(app.loc(1,app.loc_idx), app.loc(2,app.loc_idx), 'r*', 'Parent', app.ax1); % [xi, yi]
                end
                
                if app.loc_idx == 2
                    uiresume
                end
                
            elseif (app.rect_or_pixel == 2 && app.PixelRemoveButton.Value == 1) ||...
                (app.rect_or_pixel == 2 && app.PixelAddButton.Value == 1) 
                % pixel edit status
                % One point for pixel
                % keep the state as much as you want
                
                app.loc_idx = app.loc_idx + 1;
                hold(app.ax1, 'on')
                temp = app.ax1.CurrentPoint; % Returns 2x3 array of points
                app.loc(:,app.loc_idx) = [temp(1,1); temp(1,2)]; % Gets the (x,y) coordinates % any problems here?
                plot(app.loc(1,app.loc_idx), app.loc(2,app.loc_idx), 'b*', 'Parent', app.ax1);
                
            end
            
        end
        % Button pushed function: RectRemoveButton
        %function RectRemoveButtonPushed(app, event)
            
        %end

        % Button pushed function: PixelRemoveButton
        function PixelRemoveButtonValueChanged(app, event)
             value = app.PixelRemoveButton.Value;
             % app.PixelAddButton.Value = 0; % mutual exclusive to Pixel Add
             changeOnHomebrewBtn(app);
            
             if value == 1
                 app.UIFigure.Pointer = "crosshair";
                 app.loc = zeros(2,1);
                 app.loc_idx = 0;
                 app.rect_or_pixel = 2;
                 app.window_button_down = 0;
                 uiwait;
                 
                 if app.window_button_down == 3
                     app.pixel_mask_rm_count = app.pixel_mask_rm_count + 1; % before wait? No!
                     
                     hold(app.ax1, 'off');
                     pixel_mask_single = DrawPixelMask_GUI2020(app.loc, app.img_size);
                     
                     app.pixel_mask_rm = pixel_mask_single;
                     drawPixel(app);
                     
                     checkMode_updateFigures(app);
                     
                     % Remove skipped slice (if there is)
                     [skip_check, skip_match_idx] = find(app.slice_idx == app.skip_slice_idx);
                     if skip_check
                         app.skip_slice_idx(skip_match_idx) = [];
                     end
                     
                     app.SliceskipedEditField.Value = mat2str(app.skip_slice_idx);
                 end
                 % Return mouse to arrow
                 app.UIFigure.Pointer = "arrow";
             elseif value == 0
                 % Only resume when state changed
                 uiresume;
             end
        end

        % Button pushed function: ROIRemoveButton
        function ROIRemoveButtonPushed(app, event)
             % Please make sure roi is drawn
            app.window_button_down = 5; % ROI remove
            changeOnHomebrewBtn(app);
            
             if app.window_button_down == 5
                app.roi_mask_rm_count = app.roi_mask_rm_count + 1; % before wait
                
                hold(app.ax1, 'off');
                roi = drawpolygon(app.ax1);
                app.roi_mask_rm = createMask(roi); 
                
                % Based on what mode you are on
                drawROI(app);
                checkMode_updateFigures(app);

                 % Remove skipped slice (if there is)
                [skip_check, skip_match_idx] = find(app.slice_idx == app.skip_slice_idx);
                if skip_check
                    app.skip_slice_idx(skip_match_idx) = [];
                end
                
                app.SliceskippedEditField.Value = mat2str(app.skip_slice_idx);
             end
        end
        
        % Button pushed function: RectAddButton
        %function RectAddButtonPushed(app, event)
        
        %end
        
        % Button pushed function: PixelAddButton
        function PixelAddButtonValueChanged(app, event)
            value = app.PixelAddButton.Value;
            % app.PixelRemoveButton.Value = 0; % mutual exclusive to Pixel Remove
            changeOnHomebrewBtn(app);
            
            if value == 1
                app.UIFigure.Pointer = "crosshair";
                app.loc = zeros(2,1);
                app.loc_idx = 0;
                app.rect_or_pixel = 2;
                app.window_button_down = 0;
                uiwait;
                
                if app.window_button_down == 4
                    app.pixel_mask_add_count = app.pixel_mask_add_count + 1; % before wait? No!
                    
                    hold(app.ax1, 'off');
                    pixel_mask_single = DrawPixelMask_GUI2020(app.loc, app.img_size);
                    app.pixel_mask_add = pixel_mask_single;
                    
                    drawPixel(app);
                    
                    checkMode_updateFigures(app);
                    
                    % Remove skipped slice (if there is)
                    [skip_check, skip_match_idx] = find(app.slice_idx == app.skip_slice_idx);
                    if skip_check
                        app.skip_slice_idx(skip_match_idx) = [];
                    end
                    
                    app.SliceskipedEditField.Value = mat2str(app.skip_slice_idx);
                end
                % Return mouse to arrow
                app.UIFigure.Pointer = "arrow";
                
            elseif value == 0
                % Only resume when state changed
                uiresume;
            end
        end

        % Button pushed function: ROIAddButton
        function ROIAddButtonPushed(app, event)
            app.window_button_down = 6; % ROI Add
            changeOnHomebrewBtn(app);
            
            if app.window_button_down == 6
                app.roi_mask_add_count = app.roi_mask_add_count + 1; % before wait
                
                hold(app.ax1, 'off');
                roi = drawpolygon(app.ax1);
                app.roi_mask_add = createMask(roi); 
                
                
                % Based on what mode you are on
                drawROI(app);
                
                checkMode_updateFigures(app);
                
                 % Remove skipped slice (if there is)
                [skip_check, skip_match_idx] = find(app.slice_idx == app.skip_slice_idx);
                if skip_check
                    app.skip_slice_idx(skip_match_idx) = [];
                end
                
                app.SliceskippedEditField.Value = mat2str(app.skip_slice_idx);
            end
        end

        % Button pushed function: RedoButton
        function RedoButtonPushed(app, event)
            
            app.endo_mask_moded(:,:,app.slice_idx) = app.endo_mask(:,:,app.slice_idx);
            app.epi_mask_moded(:,:,app.slice_idx) = app.epi_mask(:,:,app.slice_idx);
            app.exclusion_mask_moded(:,:,app.slice_idx) = app.exclusion_mask(:,:,app.slice_idx);
            app.remote_mask_moded(:,:,app.slice_idx) = app.remote_mask(:,:,app.slice_idx);
            app.reflow_mask_moded(:,:,app.slice_idx) = app.reflow_mask(:,:,app.slice_idx);
            
            app.myo_mask(:,:,app.slice_idx) = app.epi_mask_moded(:,:,app.slice_idx) - app.endo_mask_moded(:,:,app.slice_idx);
            app.mi_mask(:,:,app.slice_idx) = zeros(app.img_size(1), app.img_size(2));
            app.mi_mask_moded(:,:,app.slice_idx) = app.mi_mask(:,:,app.slice_idx);
            app.myo_mask_moded(:,:,app.slice_idx) = app.myo_mask(:,:,app.slice_idx);
            
            update_composite_mask(app);
            update_all_axes(app);
            refreshplot(app);
            
            % Redo also clears images saved status
            app.ImageSavedStatusLamp.Color = [1 0 0];
            
            % To make buttons mutually exclusive
            app.window_button_down = 0;
            app.rect_mask_button_down = 0;
            app.rect_mask_add_button_down = 0;
            app.rect_or_pixel = 0;
            
            app.rect_mask_count = 0; % Reset rect_mask_count
            app.rect_mask_add_count = 0;
            app.pixel_mask_rm_count = 0;
            app.pixel_mask_add_count = 0;
            app.roi_mask_rm_count = 0;
            app.roi_mask_add_count = 0;
            
            app.roi_mask_rm = zeros(app.img_size(1), app.img_size(2));
            app.roi_mask_add = zeros(app.img_size(1), app.img_size(2));
            
            % Compute button initialize
            app.computebtn = 0;
            
            app.Slider.Enable = 'off';
        end

        % Button pushed function: SkipButton
        function SkipButtonPushed(app, event)
            % everything on the go now is skipped
            app.endo_mask_moded(:,:,app.slice_idx) = app.endo_mask(:,:,app.slice_idx);
            app.epi_mask_moded(:,:,app.slice_idx) = app.epi_mask(:,:,app.slice_idx);
            app.exclusion_mask_moded(:,:,app.slice_idx) = app.exclusion_mask(:,:,app.slice_idx);
            app.remote_mask_moded(:,:,app.slice_idx) = app.remote_mask(:,:,app.slice_idx);
            app.reflow_mask_moded(:,:,app.slice_idx) = app.reflow_mask(:,:,app.slice_idx);
            
            app.myo_mask(:,:,app.slice_idx) = app.epi_mask_moded(:,:,app.slice_idx) - app.endo_mask_moded(:,:,app.slice_idx);
            app.mi_mask(:,:,app.slice_idx) = zeros(app.img_size(1), app.img_size(2));
            app.mi_mask_moded(:,:,app.slice_idx) = app.mi_mask(:,:,app.slice_idx);
            app.myo_mask_moded(:,:,app.slice_idx) = app.myo_mask(:,:,app.slice_idx);
            
            app.skip_slice_idx(end+1) = app.slice_idx;
            app.skip_slice_idx = sort(app.skip_slice_idx);
            app.SliceskippedEditField.Value = mat2str(app.skip_slice_idx);
        end

        % Selection changed function: ModeButtonGroup
        function ModeButtonGroupSelectionChanged(app, event)
            selectedButton = app.ModeButtonGroup.SelectedObject;
            app.mode_btn = selectedButton.Text;
        end
        
        % Value changed function: AlgorithmsDropDown
        function AlgorithmsDropDownValueChanged(app, event)
            app.algo = app.AlgorithmsDropDown.Value;
            if strcmp(app.algo, 'Manual')
                app.Slider.Enable = 'on';
                myo_masked = app.img_3D(:,:,app.slice_idx) .* app.myo_mask(:,:,app.slice_idx);
                app.Slider.Limits = [0 max(myo_masked(:))];
            end
        end
        
        % Value changing function: Slider
        function SliderValueChanging(app, event)
            app.manual_thresh(app.slice_idx) = event.Value;
        end
        
        % Button pushed function: ComputeButton
        function ComputeButtonPushed(app, event)
            app.computebtn = 1;
            
            % Check if remote is avaialbe for display
            if any(any(app.remote_mask_moded(:,:,app.slice_idx)))
                executeAlgo(app);
            end
        end
        
         % Value changing function: Slider
        function UpperSliderValueChanging(app, event)
            changingValue = event.Value;
            app.upper_window = changingValue;
            
            update_windowing(app);
        end

        % Value changing function: LowerSlider
        function LowerSliderValueChanging(app, event)
            changingValue = event.Value;
            app.lower_window = changingValue;
            
            update_windowing(app);
        end
        
        function update_windowing(app, event)
            app.img_for_window = app.img_3D_rgb;
            app.img_for_window = (app.img_for_window - app.lower_window) / (app.upper_window - app.lower_window);
            app.img_for_window(app.img_for_window >= 1) = 1;
            app.img_for_window(app.img_for_window <= 0) = 0;
            
            %app.img_3D_rgb
            app.ax1_rgb(:,:,1) = app.rgb_scale*app.img_for_window + app.mi_mask_moded(:,:,app.slice_idx);
            app.ax1_rgb(:,:,2) = app.rgb_scale*app.img_for_window + app.remote_mask_moded(:,:,app.slice_idx);
            % edg  = edge(squeeze(app.myo_mask(:,:,app.slice_idx)),'Canny');
            app.ax1_rgb(:,:,3) = app.rgb_scale*app.img_for_window + app.myo_mask_moded(:,:,app.slice_idx);
            
            refreshplot(app);
        end

        % Private functions below
        function refreshplot(app, event)
            imagesc(app.ax1_rgb, 'Parent', app.ax1);
            set(app.ax1, 'visible', 'off');
            axis(app.ax1, 'image');
            
            % app.img_3D_rgb = adapthisteq(app.img_3D(:,:,app.slice_idx) / max(max(app.img_3D(:,:,app.slice_idx))));
            app.img_3D_rgb = app.img_3D(:,:,app.slice_idx) / max(max(app.img_3D(:,:,app.slice_idx)));
            imagesc(app.img_3D_rgb, 'Parent', app.ax2);
            set(app.ax2, 'visible', 'off');
            axis(app.ax2, 'image');
            
            imagesc(app.composite_mask(:,:,app.slice_idx), 'Parent', app.ax3);
            set(app.ax3, 'visible', 'off');
            axis(app.ax3, 'image');
        end
        
        function drawROI(app, event)
            
            switch app.mode_btn
                case {'Endo'}
                    app.first_edit_end_label = app.first_edit_end_label + 1;
                    
                    if app.first_edit_end_label == 1
                        app.endo_mask_moded(:,:,app.slice_idx) = (app.roi_mask_add + app.endo_mask(:,:,app.slice_idx) - app.roi_mask_rm) >= 1;
                    else
                        app.endo_mask_moded(:,:,app.slice_idx) = (app.roi_mask_add + app.endo_mask_moded(:,:,app.slice_idx) - app.roi_mask_rm) >= 1;
                    end
                    
                case {'Epi'}
                    app.first_edit_epi_label = app.first_edit_epi_label + 1;
                    
                    if app.first_edit_epi_label == 1
                        app.epi_mask_moded(:,:,app.slice_idx) = (app.roi_mask_add + app.epi_mask(:,:,app.slice_idx) - app.roi_mask_rm) >= 1;
                    else
                        app.epi_mask_moded(:,:,app.slice_idx) = (app.roi_mask_add + app.epi_mask_moded(:,:,app.slice_idx) - app.roi_mask_rm) >= 1;
                    end
                    
                case {'Exclusion'}
                    app.first_edit_exc_label = app.first_edit_exc_label + 1;
                    
                    if app.first_edit_exc_label == 1
                        app.exclusion_mask_moded(:,:,app.slice_idx) = (app.roi_mask_add + app.exclusion_mask(:,:,app.slice_idx) - app.roi_mask_rm) >= 1;
                    else
                        app.exclusion_mask_moded(:,:,app.slice_idx) = (app.roi_mask_add + app.exclusion_mask_moded(:,:,app.slice_idx) - app.roi_mask_rm) >= 1;
                    end
                    
                case {'Remote'}
                    app.first_edit_rem_label = app.first_edit_rem_label + 1;
                    
                    if app.first_edit_rem_label == 1
                        app.remote_mask_moded(:,:,app.slice_idx) = (app.roi_mask_add + app.remote_mask(:,:,app.slice_idx) - app.roi_mask_rm) >= 1;
                    else
                        app.remote_mask_moded(:,:,app.slice_idx) = (app.roi_mask_add + app.remote_mask_moded(:,:,app.slice_idx) - app.roi_mask_rm) >= 1;
                    end
                    
                case {'Reflow'}
                    app.first_edit_ref_label = app.first_edit_ref_label + 1;
                    
                    if app.first_edit_rem_label == 1
                        app.reflow_mask_moded(:,:,app.slice_idx) = (app.roi_mask_add + app.reflow_mask(:,:,app.slice_idx) - app.roi_mask_rm) >= 1;
                    else
                        app.reflow_mask_moded(:,:,app.slice_idx) = (app.roi_mask_add + app.reflow_mask_moded(:,:,app.slice_idx) - app.roi_mask_rm) >= 1;
                    end
                    app.reflow_mask_moded(:,:,app.slice_idx) = app.reflow_mask_moded(:,:,app.slice_idx) .* app.myo_mask(:,:,app.slice_idx);
            end
            
            % Update composite masks
            update_composite_mask(app);
        end
        
        function drawPixel(app, event)
            switch app.mode_btn
                case {'Endo'}
                    app.first_edit_end_label = app.first_edit_end_label + 1;
                    
                    if app.first_edit_end_label == 1
                        app.endo_mask_moded(:,:,app.slice_idx) = (app.pixel_mask_add + app.endo_mask(:,:,app.slice_idx) - app.pixel_mask_rm) >= 1;
                    else
                        app.endo_mask_moded(:,:,app.slice_idx) = (app.pixel_mask_add + app.endo_mask_moded(:,:,app.slice_idx) - app.pixel_mask_rm) >= 1;
                    end
                    
                case {'Epi'}
                    app.first_edit_epi_label = app.first_edit_epi_label + 1;
                    
                    if app.first_edit_epi_label == 1
                        app.epi_mask_moded(:,:,app.slice_idx) = (app.pixel_mask_add + app.epi_mask(:,:,app.slice_idx) - app.pixel_mask_rm) >= 1;
                    else
                        app.epi_mask_moded(:,:,app.slice_idx) = (app.pixel_mask_add + app.epi_mask_moded(:,:,app.slice_idx) - app.pixel_mask_rm) >= 1;
                    end
                    
                case {'Exclusion'}
                    app.first_edit_exc_label = app.first_edit_exc_label + 1;
                    
                    if app.first_edit_exc_label == 1
                        app.exclusion_mask_moded(:,:,app.slice_idx) = (app.pixel_mask_add + app.exclusion_mask(:,:,app.slice_idx) - app.pixel_mask_rm) >= 1;
                    else
                        app.exclusion_mask_moded(:,:,app.slice_idx) = (app.pixel_mask_add + app.exclusion_mask_moded(:,:,app.slice_idx) - app.pixel_mask_rm) >= 1;
                    end
                    
                case {'Remote'}
                    app.first_edit_rem_label = app.first_edit_rem_label + 1;
                    
                    if app.first_edit_rem_label == 1
                        app.remote_mask_moded(:,:,app.slice_idx) = (app.pixel_mask_add + app.remote_mask(:,:,app.slice_idx) - app.pixel_mask_rm) >= 1;
                    else
                        app.remote_mask_moded(:,:,app.slice_idx) = (app.pixel_mask_add + app.remote_mask_moded(:,:,app.slice_idx) - app.pixel_mask_rm) >= 1;
                    end
                case {'Reflow'}
                    app.first_edit_ref_label = app.first_edit_ref_label + 1;
                    
                    if app.first_edit_rem_label == 1
                        app.reflow_mask_moded(:,:,app.slice_idx) = (app.pixel_mask_add + app.reflow_mask(:,:,app.slice_idx) - app.roi_mask_rm) >= 1;
                    else
                        app.reflow_mask_moded(:,:,app.slice_idx) = (app.pixel_mask_add + app.reflow_mask_moded(:,:,app.slice_idx) - app.roi_mask_rm) >= 1;
                    end
                    app.reflow_mask_moded(:,:,app.slice_idx) = app.reflow_mask_moded(:,:,app.slice_idx) .* app.myo_mask(:,:,app.slice_idx);
            end
            
            % Update composite masks
            update_composite_mask(app);
        end
        
        function executeAlgo(app, event)
            if any(any(app.myo_mask(:,:,app.slice_idx)))
                switch app.algo
                    case {'Mean+5SD'}
                        temp_remote = app.remote_mask_moded(:,:,app.slice_idx) .* app.img_3D(:,:,app.slice_idx);
                        thresh = mean(nonzeros(temp_remote(:))) + 5*std(nonzeros(temp_remote(:)));
                        app.mi_mask(:,:,app.slice_idx) = (app.img_3D(:,:,app.slice_idx) > thresh) .* app.myo_mask(:,:,app.slice_idx);
                        app.mi_mask_moded(:,:,app.slice_idx) = app.mi_mask(:,:,app.slice_idx);
                        
                        app.ax1_rgb(:,:,1) = app.rgb_scale*app.img_for_window + app.mi_mask(:,:,app.slice_idx);
                        refreshplot(app);
                    case {'Mean+2SD'}
                        temp_remote = app.remote_mask_moded(:,:,app.slice_idx) .* app.img_3D(:,:,app.slice_idx);
                        thresh = mean(nonzeros(temp_remote(:))) - 2*std(nonzeros(temp_remote(:)));
                        app.mi_mask(:,:,app.slice_idx) = (app.img_3D(:,:,app.slice_idx) < thresh) .* app.myo_mask(:,:,app.slice_idx);
                        app.mi_mask_moded(:,:,app.slice_idx) = app.mi_mask(:,:,app.slice_idx);
                        
                        app.ax1_rgb(:,:,1) = app.rgb_scale*app.img_for_window + app.mi_mask(:,:,app.slice_idx);
                        refreshplot(app);
                    case {'Manual'}
                        app.mi_mask(:,:,app.slice_idx) = (app.img_3D(:,:,app.slice_idx) > app.manual_thresh(app.slice_idx)) .* app.myo_mask(:,:,app.slice_idx);
                        app.mi_mask_moded(:,:,app.slice_idx) = app.mi_mask(:,:,app.slice_idx);
                        app.ax1_rgb(:,:,1) = app.rgb_scale*app.img_for_window + app.mi_mask(:,:,app.slice_idx);
                        refreshplot(app);
                        
                end
            end
        end
        
        function mkdir_in_save_dir(app)
            app.blood_dir = cat(2, app.save_dir, 'BloodPool/');
            if ~exist(app.blood_dir, 'dir')
               mkdir(app.blood_dir); 
            end
            
            app.exclude_dir = cat(2, app.save_dir, 'excludeArea/');
            if ~exist(app.exclude_dir, 'dir')
               mkdir(app.exclude_dir); 
            end
            
            app.heart_dir = cat(2, app.save_dir, 'Heart/');
            if ~exist(app.heart_dir, 'dir')
               mkdir(app.heart_dir); 
            end
            
            app.myo_dir = cat(2, app.save_dir, 'Myocardium/');
            if ~exist(app.myo_dir, 'dir')
               mkdir(app.myo_dir); 
            end
            
            app.remote_dir = cat(2, app.save_dir, 'MyoReference/');
            if ~exist(app.remote_dir, 'dir')
               mkdir(app.remote_dir); 
            end
            
            app.mi_dir = cat(2, app.save_dir, app.tp, '_MI/');
            if ~exist(app.mi_dir, 'dir')
               mkdir(app.mi_dir); 
            end
            
            app.reflow_dir = cat(2, app.save_dir, 'reflowArea/');
            if ~exist(app.reflow_dir, 'dir')
               mkdir(app.reflow_dir); 
            end
            
            app.myo_pool_dir = cat(2, app.name_dir, 'MyoReservoir/');
            if ~exist(app.myo_pool_dir, 'dir')
               mkdir(app.myo_pool_dir); 
            end
        end
        
        function changeOnHomebrewBtn(app, event)
        % complete RectAdd if it was on
            if (app.window_button_down == 0 || app.window_button_down == 1) && app.rect_mask_button_down == 1
                uiresume;
                app.rect_mask_button_down = 0;
            elseif (app.window_button_down == 0 || app.window_button_down == 2) && app.rect_mask_add_button_down == 1
                uiresume;
                app.rect_mask_add_button_down = 0;
            elseif (app.window_button_down == 0 || app.window_button_down == 3) && app.PixelRemoveButton.Value == 1
                uiresume;
                app.PixelRemoveButton.Value = 0;
            elseif (app.window_button_down == 0 || app.window_button_down == 4) && app.PixelAddButton.Value == 1
                uiresume;
                app.PixelAddButton.Value = 0;
            end
        end
        
        function update_all_axes(app, event)
            app.ax1_rgb(:,:,1) = app.rgb_scale*app.img_for_window + app.mi_mask_moded(:,:,app.slice_idx);
            app.ax1_rgb(:,:,2) = app.rgb_scale*app.img_for_window + app.remote_mask_moded(:,:,app.slice_idx);
            % edg  = edge(squeeze(app.myo_mask(:,:,app.slice_idx)),'Canny');
            app.ax1_rgb(:,:,3) = app.rgb_scale*app.img_for_window + app.myo_mask_moded(:,:,app.slice_idx);

        end
        
        function update_composite_mask(app, event)
            app.composite_mask(:,:,app.slice_idx) = app.myo_mask(:,:,app.slice_idx) + ...
                app.remote_mask_moded(:,:,app.slice_idx) + (-3) * app.exclusion_mask_moded(:,:,app.slice_idx) + ...
                5 * app.reflow_mask_moded(:,:,app.slice_idx);
        end
        
        function checkMode_updateFigures(app, event)
            % Check if myocardium is available for display
            if any(any(app.endo_mask_moded(:,:,app.slice_idx))) && any(any(app.epi_mask_moded(:,:,app.slice_idx)))
                app.myo_mask(:,:,app.slice_idx) = app.epi_mask_moded(:,:,app.slice_idx) - app.endo_mask_moded(:,:,app.slice_idx);
                % edg  = edge(squeeze(app.myo_mask(:,:,app.slice_idx)),'Canny');
                app.ax1_rgb(:,:,3) = app.rgb_scale*app.img_for_window + app.myo_mask_moded(:,:,app.slice_idx);
                
                check_Manual_Compute(app);
                % save myocardium once it's generated.
                % No. 
                % save myocardium after every slice has it
            elseif any(any(app.myo_mask(:,:,app.slice_idx))) && any(any(app.endo_mask_moded(:,:,app.slice_idx)))
                app.epi_mask_moded(:,:,app.slice_idx) = imfill(app.myo_mask(:,:,app.slice_idx), 'holes');
                app.myo_mask(:,:,app.slice_idx) = app.epi_mask_moded(:,:,app.slice_idx) - app.endo_mask_moded(:,:,app.slice_idx);
                % edg  = edge(squeeze(app.myo_mask(:,:,app.slice_idx)),'Canny');
                app.ax1_rgb(:,:,3) = app.rgb_scale*app.img_for_window + app.myo_mask_moded(:,:,app.slice_idx);
                
                check_Manual_Compute(app);
            elseif any(any(app.myo_mask(:,:,app.slice_idx))) && any(any(app.epi_mask_moded(:,:,app.slice_idx)))
                app.endo_mask_moded(:,:,app.slice_idx) = imfill(app.myo_mask(:,:,app.slice_idx), 'holes') - app.myo_mask(:,:,app.slice_idx);
                app.myo_mask(:,:,app.slice_idx) = app.epi_mask_moded(:,:,app.slice_idx) - app.endo_mask_moded(:,:,app.slice_idx);
                % edg  = edge(squeeze(app.myo_mask(:,:,app.slice_idx)),'Canny');
                app.ax1_rgb(:,:,3) = app.rgb_scale*app.img_for_window + app.myo_mask_moded(:,:,app.slice_idx);
                
                check_Manual_Compute(app);
            end
            
            % Check if remote reference is available for display
            if any(any(app.remote_mask_moded(:,:,app.slice_idx)))
                app.ax1_rgb(:,:,2) = app.rgb_scale*app.img_for_window + app.remote_mask_moded(:,:,app.slice_idx);
            end
            
            % Check if exclude contour is available for display
            if any(any(app.exclusion_mask_moded(:,:,app.slice_idx)))
                app.mi_mask_moded(:,:,app.slice_idx) = app.mi_mask_moded(:,:,app.slice_idx) .* ~app.exclusion_mask_moded(:,:,app.slice_idx);
                app.ax1_rgb(:,:,1) = app.rgb_scale*app.img_for_window + app.mi_mask_moded(:,:,app.slice_idx);
            end
            
             % Check if reflow contour is available for display
            if any(any(app.mi_mask_moded(:,:,app.slice_idx)))
                app.mi_mask_moded(:,:,app.slice_idx) = (app.mi_mask_moded(:,:,app.slice_idx) + app.reflow_mask_moded(:,:,app.slice_idx)) >= 1;                
                app.ax1_rgb(:,:,1) = app.rgb_scale*app.img_for_window + app.mi_mask_moded(:,:,app.slice_idx);
            end
            
            update_composite_mask(app);
            refreshplot(app);
        end
        
        function check_Manual_Compute(app, event)
            myo_masked = app.img_3D(:,:,app.slice_idx) .* app.myo_mask(:,:,app.slice_idx);
            if strcmp(app.algo, 'Manual') && any(myo_masked(:))
                app.Slider.Enable = 'on';
                app.Slider.Limits = [0 max(myo_masked(:))];
            end
        end
        
        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {735, 735};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {270, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end
            
    % Try what was supposed in addpath('functions/')
    end


    % Component initialization
    methods (Access = private)
        
        % Create UIFigure and components
        function createComponents(app)
            
            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 1224 735];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);
            app.UIFigure.WindowButtonDownFcn = createCallbackFcn(app, @UIFigureWindowButtonDown, true);
            
            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {270, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'off';
            
            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;
            app.LeftPanel.Scrollable = 'off';
            
            % Create ImageModalityDropDownLabel
            app.ImageModalityDropDownLabel = uilabel(app.LeftPanel);
            app.ImageModalityDropDownLabel.HorizontalAlignment = 'right';
            app.ImageModalityDropDownLabel.FontSize = 24;
            app.ImageModalityDropDownLabel.Position = [7 700 171 32];
            app.ImageModalityDropDownLabel.Text = 'Name List';
            
            % Create ImageModalityDropDown
            app.ImageModalityDropDown = uidropdown(app.LeftPanel);
            app.name_glob = glob(cat(2, 'PatientFatAnalysis/', 'ContourData', '/*/'));
            Names = ExtractNames(app.name_glob);
            dropdown_list = ['Please Select'; Names];
            app.ImageModalityDropDown.Items = dropdown_list;
            app.ImageModalityDropDown.ValueChangedFcn = createCallbackFcn(app, @ImageModalityDropDownValueChanged, true);
            app.ImageModalityDropDown.FontSize = 24;
            app.ImageModalityDropDown.Position = [121 662 143 31];
            
            % Create NameListDropDownLabel
            app.NameListDropDownLabel = uilabel(app.LeftPanel);
            app.NameListDropDownLabel.HorizontalAlignment = 'right';
            app.NameListDropDownLabel.FontSize = 24;
            app.NameListDropDownLabel.Position = [8 599 114 32];
            app.NameListDropDownLabel.Text = 'Time Point';
            
            % Create NameListDropDown
            app.NameListDropDown = uidropdown(app.LeftPanel);
            app.tp_glob = glob(cat(2, 'PatientFatAnalysis/', 'ContourData/', '484060000012/*/'));
            tps = ExtractNames(app.tp_glob);
            app.NameListDropDown.Items = ['Please Select'; tps];
            app.NameListDropDown.ValueChangedFcn = createCallbackFcn(app, @NameListDropDownValueChanged, true);
            app.NameListDropDown.FontSize = 24;
            app.NameListDropDown.Position = [136 600 128 31];
            
            % Create LoadFileButton
            app.LoadFileButton = uibutton(app.LeftPanel, 'push');
            app.LoadFileButton.ButtonPushedFcn = createCallbackFcn(app, @LoadFileButtonPushed, true);
            app.LoadFileButton.FontSize = 24;
            app.LoadFileButton.Position = [7 524 113 39];
            app.LoadFileButton.Text = 'Load File';
            
            % Create SaveFileButton
            app.SaveFileButton = uibutton(app.LeftPanel, 'push');
            app.SaveFileButton.ButtonPushedFcn = createCallbackFcn(app, @SaveFileButtonPushed, true);
            app.SaveFileButton.FontSize = 24;
            app.SaveFileButton.Position = [136 524 113 39];
            app.SaveFileButton.Text = 'Save File';
            
            % Create LoadCorrectedMaskButton
            app.LoadCorrectedMaskButton = uibutton(app.LeftPanel, 'push');
            app.LoadCorrectedMaskButton.ButtonPushedFcn = createCallbackFcn(app, @LoadCorrectedMaskButtonPushed, true);
            app.LoadCorrectedMaskButton.FontSize = 24;
            app.LoadCorrectedMaskButton.Position = [8 464 250 39];
            app.LoadCorrectedMaskButton.Text = 'Load Myo Mask';
            
            % Create SliceEditFieldLabel
            app.SliceEditFieldLabel = uilabel(app.LeftPanel);
            app.SliceEditFieldLabel.HorizontalAlignment = 'right';
            app.SliceEditFieldLabel.FontSize = 24;
            app.SliceEditFieldLabel.Position = [53 167 58 32];
            app.SliceEditFieldLabel.Text = 'Slice';
            
            % Create SliceEditField
            app.SliceEditField = uieditfield(app.LeftPanel, 'text');
            app.SliceEditField.FontSize = 24;
            app.SliceEditField.Position = [126 165 124 34];
            
            % Create SliceskippedEditFieldLabel
            app.SliceskippedEditFieldLabel = uilabel(app.LeftPanel);
            app.SliceskippedEditFieldLabel.HorizontalAlignment = 'right';
            app.SliceskippedEditFieldLabel.FontSize = 24;
            app.SliceskippedEditFieldLabel.Position = [5 110 149 32];
            app.SliceskippedEditFieldLabel.Text = 'Slice skipped';
            
            % Create SliceskippedEditField
            app.SliceskippedEditField = uieditfield(app.LeftPanel, 'text');
            app.SliceskippedEditField.FontSize = 24;
            app.SliceskippedEditField.Position = [169 108 81 33.8000011444092];
            
            % Create ImageSavedStatusLampLabel
            app.ImageSavedStatusLampLabel = uilabel(app.LeftPanel);
            app.ImageSavedStatusLampLabel.HorizontalAlignment = 'right';
            app.ImageSavedStatusLampLabel.FontSize = 24;
            app.ImageSavedStatusLampLabel.Position = [4 55 222 32];
            app.ImageSavedStatusLampLabel.Text = 'Image Saved Status';
            
            % Create ImageSavedStatusLamp
            app.ImageSavedStatusLamp = uilamp(app.LeftPanel);
            app.ImageSavedStatusLamp.Position = [237 59 26 26];
            
            % Create WindowingPanel
            app.WindowingPanel = uipanel(app.LeftPanel);
            app.WindowingPanel.Title = 'Windowing';
            app.WindowingPanel.FontSize = 20;
            app.WindowingPanel.Position = [5 236 259 211];

            % Create UpperSliderLabel
            app.UpperSliderLabel = uilabel(app.WindowingPanel);
            app.UpperSliderLabel.HorizontalAlignment = 'right';
            app.UpperSliderLabel.FontSize = 18;
            app.UpperSliderLabel.Position = [11 135 55 23];
            app.UpperSliderLabel.Text = 'Upper';

            % Create UpperSlider
            app.UpperSlider = uislider(app.WindowingPanel);
            app.UpperSlider.ValueChangingFcn = createCallbackFcn(app, @UpperSliderValueChanging, true);
            app.UpperSlider.FontSize = 18;
            app.UpperSlider.Position = [86 145 150 3];

            % Create LowerSliderLabel
            app.LowerSliderLabel = uilabel(app.WindowingPanel);
            app.LowerSliderLabel.HorizontalAlignment = 'right';
            app.LowerSliderLabel.FontSize = 18;
            app.LowerSliderLabel.Position = [10 50 55 23];
            app.LowerSliderLabel.Text = 'Lower';

            % Create LowerSlider
            app.LowerSlider = uislider(app.WindowingPanel);
            app.LowerSlider.ValueChangingFcn = createCallbackFcn(app, @LowerSliderValueChanging, true);
            app.LowerSlider.FontSize = 18;
            app.LowerSlider.Position = [86 60 150 3];
            
            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;
            app.RightPanel.Scrollable = 'off';
            
            % Create OrigImagePanel
            app.OrigImagePanel = uipanel(app.RightPanel);
            app.OrigImagePanel.Title = 'Orig Image';
            app.OrigImagePanel.FontSize = 20;
            app.OrigImagePanel.Position = [4 236 488.99 412];
            
            % Create EditPanel
            app.EditPanel = uipanel(app.RightPanel);
            app.EditPanel.Title = 'Edit Panel';
            app.EditPanel.FontSize = 20;
            app.EditPanel.Position = [716 236 238 496];
            
            % Create RectRemoveButton
            %app.RectRemoveButton = uibutton(app.EditPanel, 'push');
            %app.RectRemoveButton.ButtonPushedFcn = createCallbackFcn(app, @RectRemoveButtonPushed, true);
            %app.RectRemoveButton.FontSize = 22;
            %app.RectRemoveButton.Position = [34 405 171 39];
            %app.RectRemoveButton.Text = 'Rect Remove';
            
            % Create PixelRemoveButton
            app.PixelRemoveButton = uibutton(app.EditPanel, 'state');
            app.PixelRemoveButton.ValueChangedFcn = createCallbackFcn(app, @PixelRemoveButtonValueChanged, true);
            app.PixelRemoveButton.FontSize = 22;
            app.PixelRemoveButton.Position = [34 404 170 39];
            app.PixelRemoveButton.Text = 'Pixel Remove';
            
            % Create ROIRemoveButton
            app.ROIRemoveButton = uibutton(app.EditPanel, 'push');
            app.ROIRemoveButton.ButtonPushedFcn = createCallbackFcn(app, @ROIRemoveButtonPushed, true);
            app.ROIRemoveButton.FontSize = 22;
            app.ROIRemoveButton.Position = [35 355 170 39];
            app.ROIRemoveButton.Text = 'ROI Remove';
            
            % Create RectAddButton
            %app.RectAddButton = uibutton(app.EditPanel, 'push');
            %app.RectAddButton.ButtonPushedFcn = createCallbackFcn(app, @RectAddButtonPushed, true);
            %app.RectAddButton.FontSize = 22;
            %app.RectAddButton.Position = [35 246 170 39];
            %app.RectAddButton.Text = 'Rect Add';
            
            % Create PixelAddButton
            app.PixelAddButton = uibutton(app.EditPanel, 'state');
            app.PixelAddButton.ValueChangedFcn = createCallbackFcn(app, @PixelAddButtonValueChanged, true);
            app.PixelAddButton.FontSize = 22;
            app.PixelAddButton.Position = [35 288 170 39];
            app.PixelAddButton.Text = 'Pixel Add';
            
            % Create ROIAddButton
            app.ROIAddButton = uibutton(app.EditPanel, 'push');
            app.ROIAddButton.ButtonPushedFcn = createCallbackFcn(app, @ROIAddButtonPushed, true);
            app.ROIAddButton.FontSize = 22;
            app.ROIAddButton.Position = [35 239 170 39];
            app.ROIAddButton.Text = 'ROI Add';
            
            % Create RedoButton
            app.RedoButton = uibutton(app.EditPanel, 'push');
            app.RedoButton.ButtonPushedFcn = createCallbackFcn(app, @RedoButtonPushed, true);
            app.RedoButton.FontSize = 22;
            app.RedoButton.Position = [35 171 170 39];
            app.RedoButton.Text = 'Redo';
            
            % Create SkipButton
            app.SkipButton = uibutton(app.EditPanel, 'push');
            app.SkipButton.ButtonPushedFcn = createCallbackFcn(app, @SkipButtonPushed, true);
            app.SkipButton.FontSize = 22;
            app.SkipButton.Position = [35 121 170 39];
            app.SkipButton.Text = 'Skip';
            
            % Create ModeButtonGroup
            app.ModeButtonGroup = uibuttongroup(app.RightPanel);
            app.ModeButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @ModeButtonGroupSelectionChanged, true);
            app.ModeButtonGroup.Title = 'Mode';
            app.ModeButtonGroup.FontSize = 20;
            app.ModeButtonGroup.Position = [4 647 488.99 84.98];
            
            % Create EndoButton
            app.EndoButton = uitogglebutton(app.ModeButtonGroup);
            app.EndoButton.Text = 'Endo';
            app.EndoButton.FontSize = 20;
            app.EndoButton.Position = [10 13 64 32];
            app.EndoButton.Value = true;
            
            % Create EpiButton
            app.EpiButton = uitogglebutton(app.ModeButtonGroup);
            app.EpiButton.Text = 'Epi';
            app.EpiButton.FontSize = 20;
            app.EpiButton.Position = [85 13 64 32];
            
            % Create ExclusionButton
            app.ExclusionButton = uitogglebutton(app.ModeButtonGroup);
            app.ExclusionButton.Text = 'Exclusion';
            app.ExclusionButton.FontSize = 20;
            app.ExclusionButton.Position = [160 13 101 32];
            
            % Create RemoteButton
            app.RemoteButton = uitogglebutton(app.ModeButtonGroup);
            app.RemoteButton.Text = 'Remote';
            app.RemoteButton.FontSize = 20;
            app.RemoteButton.Position = [274 13 98 32];
            
            % Create ReflowButton
            app.ReflowButton = uitogglebutton(app.ModeButtonGroup);
            app.ReflowButton.Text = 'Reflow';
            app.ReflowButton.FontSize = 20;
            app.ReflowButton.Position = [384 13 97 32];
            
            % Create Panel_2
            app.Panel_2 = uipanel(app.RightPanel);
            app.Panel_2.Position = [716 6 237 230];
            
            % Create PreviousButton_2
            app.PreviousButton_2 = uibutton(app.EditPanel, 'push');
            app.PreviousButton_2.ButtonPushedFcn = createCallbackFcn(app, @PreviousButton_2Pushed, true);
            app.PreviousButton_2.FontSize = 20;
            app.PreviousButton_2.Position = [11 41 101 39];
            app.PreviousButton_2.Text = 'Previous';
            
            % Create NextButton_2
            app.NextButton_2 = uibutton(app.EditPanel, 'push');
            app.NextButton_2.ButtonPushedFcn = createCallbackFcn(app, @NextButton_2Pushed, true);
            app.NextButton_2.FontSize = 20;
            app.NextButton_2.Position = [128 41 101 39];
            app.NextButton_2.Text = 'Next';
            
            % Create AlgorithmsDropDownLabel
            app.AlgorithmsDropDownLabel = uilabel(app.Panel_2);
            app.AlgorithmsDropDownLabel.HorizontalAlignment = 'right';
            app.AlgorithmsDropDownLabel.FontSize = 20;
            app.AlgorithmsDropDownLabel.Position = [11 192 101 26];
            app.AlgorithmsDropDownLabel.Text = 'Algorithms';
            
            % Create AlgorithmsDropDown
            app.AlgorithmsDropDown = uidropdown(app.Panel_2);
            app.AlgorithmsDropDown.Items = ['Please Select One Algorithm', app.algo_label];
            app.AlgorithmsDropDown.ValueChangedFcn = createCallbackFcn(app, @AlgorithmsDropDownValueChanged, true);
            app.AlgorithmsDropDown.FontSize = 20;
            app.AlgorithmsDropDown.Position = [34 153 184 25];
            
            % Create Slider
            app.Slider = uislider(app.Panel_2);
            app.Slider.ValueChangedFcn = createCallbackFcn(app, @refreshplot, true);
            app.Slider.ValueChangingFcn = createCallbackFcn(app, @SliderValueChanging, true);
            app.Slider.FontSize = 20;
            app.Slider.Position = [51 118 150 3];
            
            % Create ComputeButton
            app.ComputeButton = uibutton(app.Panel_2, 'push');
            app.ComputeButton.ButtonPushedFcn = createCallbackFcn(app, @ComputeButtonPushed, true);
            app.ComputeButton.FontSize = 20;
            app.ComputeButton.Position = [66 18 100 32];
            app.ComputeButton.Text = 'Compute';
            
            % Create ImageOverviewPanel
            %app.ImageOverviewPanelGroup = uitabgroup(fig,'Position',[20 20 196 145]);
            %app.ImageOverviewPanel = uitab(tg,'Title','Member Information');
            app.ImageOverviewPanel = uipanel(app.RightPanel);
            app.ImageOverviewPanel.Title = 'Image Overview';
            app.ImageOverviewPanel.Scrollable = 'on';
            app.ImageOverviewPanel.FontSize = 20;
            app.ImageOverviewPanel.Position = [492 6 225 726];
            
             % Create TabGroup
            app.TabGroup = uitabgroup(app.ImageOverviewPanel);
            app.TabGroup.Position = [1 -1 224 695];

            % Create Tab
            app.Tab = uitab(app.TabGroup);
            app.Tab.Title = 'Slice 1-4';
            %app.Tab.Scrollable = 'on';
            
            app.Tab2 = uitab(app.TabGroup);
            app.Tab2.Title = 'Slice 5-8';
            
            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end
    
    % App creation and deletion
    methods (Access = public)
        
        % Construct app
        function app = GUI2022
            if (~isdeployed)
                app.funcPath = addpath('functions/'); % Addpath to functions
            end
            
            % Create UIFigure and components
            createComponents(app)
            
            % Register the app with App Designer
            registerApp(app, app.UIFigure)
            
            % Execute the startup function
            runStartupFcn(app, @startupFcn)
            
            if nargout == 0
                clear app
            end
        end
        
        % Code that executes before app deletion
        function delete(app)
            
            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end