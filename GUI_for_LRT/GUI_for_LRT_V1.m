classdef GUI_for_LRT_V1 < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                matlab.ui.Figure
        UIAxes                  matlab.ui.control.UIAxes
        ControlPanel            matlab.ui.container.Panel
        SliceSliderLabel        matlab.ui.control.Label
        SliceSlider             matlab.ui.control.Slider
        T1SliderLabel           matlab.ui.control.Label
        T1Slider                matlab.ui.control.Slider
        T2starSliderLabel       matlab.ui.control.Label
        T2starSlider            matlab.ui.control.Slider
        RespiratorySliderLabel  matlab.ui.control.Label
        RespiratorySlider       matlab.ui.control.Slider
        CardiacSliderLabel      matlab.ui.control.Label
        CardiacSlider           matlab.ui.control.Slider
        LoadDataButton          matlab.ui.control.Button
    end
    
    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        Phi
        params
        Nx
        Ny
        Nz
        L
        U
        Gr
        
        ax1
        vec
        dispim
        img_label % 0 - default; 1 - slice; 2 - T1; 4 - resp; 3 - cardiac; 5 - T2star;
        img_label_old
        
        temp1
        cw
    end
    
    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: LoadDataButton
        function LoadDataButtonPushed(app, event)
            
            [fid_file, fid_path] = uigetfile('*.mat');
            df = load(strcat(fid_path, fid_file));
            app.L = df.L;
            app.Gr = df.Gr;
            app.U = df.U;
            app.Ny = df.Ny;
            app.Nx = df.Nx;
            app.Nz = df.Nz;
            app.params = df.params;
            
            app.Phi = reshape(df.Phi, [app.L, df.sizes(2:end)]);
            app.vec = df.vec;
            app.dispim = df.dispim;
            %slider_dim = length(df.sizes)-1;
            
            app.SliceSlider.Limits = [1, app.Nz];
            app.T1Slider.Limits = [1, df.sizes(2)];
            app.CardiacSlider.Limits = [1, df.sizes(3)];
            app.RespiratorySlider.Limits = [1, df.sizes(4)];
            app.T2starSlider.Limits = [1 df.sizes(5)];
            
            app.SliceSlider.Value = 4;
            app.T1Slider.Value = df.sizes(2);
            app.RespiratorySlider.Value = 1;
            app.CardiacSlider.Value = 1;
            app.T2starSlider.Value = 1;
            
            app.SliceSlider.MajorTicks = [1:2:app.Nz];
            app.T1Slider.MajorTicks = [1:20:df.sizes(2)];
            app.CardiacSlider.MajorTicks = [1:1:df.sizes(3)];
            app.RespiratorySlider.MajorTicks = [1:1:df.sizes(4)];
            app.T2starSlider.MajorTicks = [1:1:df.sizes(5)];
            
            %app.T1Slider.MajorTickLabels = {};
            
            app.SliceSlider.MinorTicks = [];
            app.T1Slider.MinorTicks = [];
            app.RespiratorySlider.MinorTicks = [];
            app.CardiacSlider.MinorTicks = [];
            app.T2starSlider.MinorTicks = [];
            
            app.img_label = 0;
            app.img_label_old = 0;
            app.ax1 = uiaxes(app.UIFigure,'Position', [1 1 441 480]);
            
            update_axes(app);
            
        end
        
        % Value changed function: SliceSlider
        function SliceSliderValueChanged(app, event)
            app.img_label = 1;
            value = app.SliceSlider.Value;
            slice_range = [1:1:app.Nz];
            % determine which discrete option the current value is closest to.
            [~, minIdx] = min(abs(value - event.Source.MajorTicks(:)));
            % move the slider to that option
            event.Source.Value = event.Source.MajorTicks(minIdx);
            % Override the selected value if you plan on using it within this function
            app.SliceSlider.Value = event.Source.MajorTicks(minIdx);
            
            update_axes(app);
            app.img_label_old = 1;
        end

        % Value changed function: T1Slider
        function T1SliderValueChanged(app, event)
            app.img_label = 2;
            value = app.T1Slider.Value;
            % t1_range = [1:1:app.];
            % determine which discrete option the current value is closest to.
            [~, minIdx] = min(abs(value - event.Source.MajorTicks(:)));
            % move the slider to that option
            event.Source.Value = event.Source.MajorTicks(minIdx);
            % Override the selected value if you plan on using it within this function
            app.T1Slider.Value = event.Source.MajorTicks(minIdx);
            
            update_axes(app);
            app.img_label_old = 2;
        end

        % Value changed function: T2starSlider
        function T2starSliderValueChanged(app, event)
            app.img_label = 5;
            value = app.T2starSlider.Value;
             % determine which discrete option the current value is closest to.
            [~, minIdx] = min(abs(value - event.Source.MajorTicks(:)));
            % move the slider to that option
            event.Source.Value = event.Source.MajorTicks(minIdx);
            % Override the selected value if you plan on using it within this function
            app.T2starSlider.Value = event.Source.MajorTicks(minIdx);
            
            update_axes(app);
            app.img_label_old = 5;
        end

        % Value changed function: RespiratorySlider
        function RespiratorySliderValueChanged(app, event)
            app.img_label = 4;
            value = app.RespiratorySlider.Value;
            
             % determine which discrete option the current value is closest to.
            [~, minIdx] = min(abs(value - event.Source.MajorTicks(:)));
            % move the slider to that option
            event.Source.Value = event.Source.MajorTicks(minIdx);
            % Override the selected value if you plan on using it within this function
            app.RespiratorySlider.Value = event.Source.MajorTicks(minIdx);
            
            update_axes(app);
            app.img_label_old = 4;
        end

        % Value changed function: CardiacSlider
        function CardiacSliderValueChanged(app, event)
            app.img_label = 3;
            value = app.CardiacSlider.Value;
            
             % determine which discrete option the current value is closest to.
            [~, minIdx] = min(abs(value - event.Source.MajorTicks(:)));
            % move the slider to that option
            event.Source.Value = event.Source.MajorTicks(minIdx);
            % Override the selected value if you plan on using it within this function
            app.CardiacSlider.Value = event.Source.MajorTicks(minIdx);
            
            update_axes(app);
            app.img_label_old = 3;
        end
        
        function update_axes(app, event)
            img_label_diff = abs(app.img_label - app.img_label_old);
            
            if app.img_label == 1 && img_label_diff == 0
                app.dispim = @(x)fftshift(x(:,:,app.SliceSlider.Value,:),1);
            end
            
            if (app.img_label == 2 && img_label_diff ~= 0) || (app.img_label == 5 && img_label_diff ~= 0)
                app.temp1 = app.Gr\reshape(app.Phi(:,:,app.CardiacSlider.Value,app.RespiratorySlider.Value,:), app.L, []);
                app.temp1 = reshape(reshape(app.dispim(reshape(app.U,app.Ny,app.Nx,app.Nz,[])),[],app.L)*app.temp1, app.Ny, app.Nx, [], app.params.NEco);
                app.cw = 0.5*max(app.vec(abs(app.temp1)));
            elseif app.img_label == 3 && img_label_diff ~= 0
                app.temp1 = app.Gr\reshape(app.Phi(:,app.T1Slider.Value,:,app.RespiratorySlider.Value,:), app.L, []);
                app.temp1 = reshape(reshape(app.dispim(reshape(app.U,app.Ny,app.Nx,app.Nz,[])),[],app.L)*app.temp1, app.Ny, app.Nx, [], app.params.NEco);
                app.cw = 0.5*max(app.vec(abs(app.temp1)));
            elseif app.img_label == 4 && img_label_diff ~= 0
                app.temp1 = app.Gr\reshape(app.Phi(:,app.T1Slider.Value,app.CardiacSlider.Value,:,:), app.L, []);
                app.temp1 = reshape(reshape(app.dispim(reshape(app.U,app.Ny,app.Nx,app.Nz,[])),[],app.L)*app.temp1, app.Ny, app.Nx, [], app.params.NEco);
                app.cw = 0.5*max(app.vec(abs(app.temp1)));
            elseif app.img_label == 0
                app.temp1 = app.Gr\reshape(app.Phi(:,end,1,1,:), app.L, []);
                app.temp1 = reshape(reshape(app.dispim(reshape(app.U,app.Ny,app.Nx,app.Nz,[])),[],app.L)*app.temp1, app.Ny, app.Nx, [], app.params.NEco);
                app.cw = 0.5*max(app.vec(abs(app.temp1)));
            elseif app.img_label == 1 
                app.temp1 = app.Gr\reshape(app.Phi(:,app.T1Slider.Value,app.CardiacSlider.Value,app.RespiratorySlider.Value,:), app.L, []);
                app.temp1 = reshape(reshape(app.dispim(reshape(app.U,app.Ny,app.Nx,app.Nz,[])),[],app.L)*app.temp1, app.Ny, app.Nx, [], app.params.NEco);
                app.cw = 0.5*max(app.vec(abs(app.temp1)));
            end
            
            
            if app.img_label == 2 || app.img_label == 5
                imagesc(abs(app.temp1(:,:,app.T1Slider.Value, app.T2starSlider.Value)/app.cw), 'Parent', app.ax1);
            elseif app.img_label == 3
                imagesc(abs(app.temp1(:,:,app.CardiacSlider.Value, app.T2starSlider.Value)/app.cw), 'Parent', app.ax1);
            elseif app.img_label == 4
                imagesc(abs(app.temp1(:,:,app.RespiratorySlider.Value, app.T2starSlider.Value)/app.cw), 'Parent', app.ax1);
            elseif app.img_label == 0 || app.img_label == 1
                imagesc(abs(app.temp1(:,:,:,app.T2starSlider.Value)/app.cw), 'Parent', app.ax1);
            end
            
            set(app.ax1, 'visible', 'off');
            axis(app.ax1, 'image');
            colormap(app.ax1, gray);
        end
        
%         function dispim(app, event)
%             
%         end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 686 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            app.UIAxes.PlotBoxAspectRatio = [1 1.08163265306122 1];
            app.UIAxes.Position = [1 1 465 480];

            % Create ControlPanel
            app.ControlPanel = uipanel(app.UIFigure);
            app.ControlPanel.Title = 'Control Panel';
            app.ControlPanel.Position = [465 1 222 480];

            % Create SliceSliderLabel
            app.SliceSliderLabel = uilabel(app.ControlPanel);
            app.SliceSliderLabel.HorizontalAlignment = 'right';
            app.SliceSliderLabel.Position = [41 428 31 22];
            app.SliceSliderLabel.Text = 'Slice';

            % Create SliceSlider
            app.SliceSlider = uislider(app.ControlPanel);
            app.SliceSlider.ValueChangedFcn = createCallbackFcn(app, @SliceSliderValueChanged, true);
            app.SliceSlider.Position = [90 438 117 3];

            % Create T1SliderLabel
            app.T1SliderLabel = uilabel(app.ControlPanel);
            app.T1SliderLabel.HorizontalAlignment = 'right';
            app.T1SliderLabel.Position = [47 363 25 22];
            app.T1SliderLabel.Text = 'T1';

            % Create T1Slider
            app.T1Slider = uislider(app.ControlPanel);
            app.T1Slider.ValueChangedFcn = createCallbackFcn(app, @T1SliderValueChanged, true);
            app.T1Slider.Position = [90 373 117 3];

            % Create T2starSliderLabel
            app.T2starSliderLabel = uilabel(app.ControlPanel);
            app.T2starSliderLabel.HorizontalAlignment = 'right';
            app.T2starSliderLabel.Position = [33 296 39 22];
            app.T2starSliderLabel.Text = 'T2star';

            % Create T2starSlider
            app.T2starSlider = uislider(app.ControlPanel);
            app.T2starSlider.ValueChangedFcn = createCallbackFcn(app, @T2starSliderValueChanged, true);
            app.T2starSlider.Position = [90 306 117 3];

            % Create RespiratorySliderLabel
            app.RespiratorySliderLabel = uilabel(app.ControlPanel);
            app.RespiratorySliderLabel.HorizontalAlignment = 'right';
            app.RespiratorySliderLabel.Position = [3 224 67 22];
            app.RespiratorySliderLabel.Text = 'Respiratory';

            % Create RespiratorySlider
            app.RespiratorySlider = uislider(app.ControlPanel);
            app.RespiratorySlider.ValueChangedFcn = createCallbackFcn(app, @RespiratorySliderValueChanged, true);
            app.RespiratorySlider.Position = [90 234 117 3];

            % Create cardiacSliderLabel
            app.CardiacSliderLabel = uilabel(app.ControlPanel);
            app.CardiacSliderLabel.HorizontalAlignment = 'right';
            app.CardiacSliderLabel.Position = [27 153 45 22];
            app.CardiacSliderLabel.Text = 'Cardiac';

            % Create cardiacSlider
            app.CardiacSlider = uislider(app.ControlPanel);
            app.CardiacSlider.ValueChangedFcn = createCallbackFcn(app, @CardiacSliderValueChanged, true);
            app.CardiacSlider.Position = [90 163 117 3];

            % Create LoadDataButton
            app.LoadDataButton = uibutton(app.ControlPanel, 'push');
            app.LoadDataButton.ButtonPushedFcn = createCallbackFcn(app, @LoadDataButtonPushed, true);
            app.LoadDataButton.Position = [8 84 100 22];
            app.LoadDataButton.Text = 'Load Data';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = GUI_for_LRT_V1

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

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