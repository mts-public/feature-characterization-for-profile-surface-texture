classdef GUI < matlab.apps.AppBase
% GUI for the feature characterization of profiles according to
% ISO 21920-2. The computation uses the functions of the folders
% "featurecharacterization2d" and "softgauge"; only the plot of the motifs
% into the app axes is part of the app.
% Profiles can be loaded from softgauge files (*.smd, step size dx from
% the file) or *.mat files with the variable z (in µm) and optionally dx
% (in mm, default 0.5e-3 mm). The mean of z is subtracted.

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        ConventionISO21920Panel        matlab.ui.container.Panel
        mLabel                         matlab.ui.control.Label
        Label_2                        matlab.ui.control.Label
        Label                          matlab.ui.control.Label
        FeatureCharacterizationEditField  matlab.ui.control.EditField
        FeatureparameterPanel          matlab.ui.container.Panel
        CheckBox_2                     matlab.ui.control.CheckBox
        AttributetypeDropDown          matlab.ui.control.DropDown
        AttributetypeDropDownLabel     matlab.ui.control.Label
        StatisticsvalueEditField       matlab.ui.control.NumericEditField
        StatisticsvalueEditFieldLabel  matlab.ui.control.Label
        NIsigEditField                 matlab.ui.control.NumericEditField
        NIsigEditFieldLabel            matlab.ui.control.Label
        AttributestatisticsDropDown    matlab.ui.control.DropDown
        AttributestatisticsDropDownLabel  matlab.ui.control.Label
        SignificantfeaturesDropDown    matlab.ui.control.DropDown
        SignificantfeaturesDropDownLabel  matlab.ui.control.Label
        LoadprofileButton              matlab.ui.control.Button
        ThresholdEditField             matlab.ui.control.NumericEditField
        ThresholdEditFieldLabel        matlab.ui.control.Label
        FeaturetypeDropDown            matlab.ui.control.DropDown
        FeaturetypeDropDownLabel       matlab.ui.control.Label
        PruningTypeDropDown            matlab.ui.control.DropDown
        PruningTypeDropDownLabel       matlab.ui.control.Label
        WatershedsegmentationPanel     matlab.ui.container.Panel
        CheckBox                       matlab.ui.control.CheckBox
        UIAxes                         matlab.ui.control.UIAxes
    end

    properties (Access = private)
        dx  % step size in x-direction in mm
        z   % vertical profile values in µm
    end

    methods (Access = private)
        function update_FC(app)
            FT = app.FeaturetypeDropDown.Value;
            PT = app.PruningTypeDropDown.Value;
            TH_perc = app.CheckBox.Value;
            if PT == "None"
                TH = "";
            elseif TH_perc == 0
                TH = " " + num2str(app.ThresholdEditField.Value);
            else
                TH = " " + num2str(app.ThresholdEditField.Value) + " %";
            end
            Fsig = app.SignificantfeaturesDropDown.Value;
            NIsig_perc = app.CheckBox_2.Value;
            if Fsig == "All"
                NIsig = "";
            elseif NIsig_perc == 0
                NIsig = " " + num2str(app.NIsigEditField.Value);
            else
                NIsig = " " + num2str(app.NIsigEditField.Value) + " %";
            end
            AT = app.AttributetypeDropDown.Value;
            Astats = app.AttributestatisticsDropDown.Value;
            if Astats == "Perc" || Astats == "Hist"
                vstats = " " + num2str(app.StatisticsvalueEditField.Value);
            else
                vstats = "";
            end
            app.FeatureCharacterizationEditField.Value = "FC; " + FT + "; " + ...
                PT + TH + "; " + Fsig + NIsig + "; " + AT + "; " + Astats + vstats;
        end

        function update_plot(app)
            cla(app.UIAxes)
            % no profile loaded yet
            if isempty(app.z)
                return
            end
            FC_complete = convertCharsToStrings(app.FeatureCharacterizationEditField.Value);
            FC = strtrim(split(FC_complete, ";"));
            [xFC, M, META] = feature_characterization(app.z, app.dx, FC(2), ...
                FC(3), FC(4), FC(5), FC(6));
            hold(app.UIAxes, 'on')
            plot_motifs(app, app.UIAxes, app.z, app.dx, M, META.Fsig, META.NIsig)
            if isnumeric(xFC)
                app.Label.Text = num2str(xFC, 2);
            else
                % "Hist": histogram in a separate figure
                app.Label.Text = "-";
            end
            app.mLabel.Text = parameter_unit(META.AT, META.Astats);
        end

        %% plot of the motifs into the app axes
        function plot_motifs(app, ax, z, dx, M, Fsig, NIsig)
            % INPUTS:
            %   ax    - axes to plot into
            %   z     - vertical profile values in µm
            %   dx    - step size in x-direction in mm
            %   M     - motif array
            %   Fsig  - significant features {"All", "Open", "Closed", "Top", "Bot"}
            %   NIsig - nesting index for significant features
            z = z(:);
            x = (0:length(z) - 1)'*dx;
            % plot settings
            linewidth = 0.85;
            % transparency-value for not significant motifs
            alpha = 0.5;
            % color of motifs
            col = [0 0.4470 0.7410];
            xlabel(ax, 'profile length / mm')
            ylabel(ax, 'profile height / µm')
            ax.XLim = [0 max(x)];
            ax.YLim = 1.1*[min(z) max(z)];
            ax.XGrid = 'on';
            ax.YGrid = 'on';
            ax.Box = 'on';
            % fill motifs and motif "frames"
            for k = 1:length(M)
                Mk = M(k);
                p = patch_featureelement(app, ax, z, dx, Mk, col);
                pm = plot(ax, ([Mk.ilp Mk.ilp Mk.ihp Mk.ihp] - 1)*dx, ...
                    [z(ceil(Mk.ilp)) z(ceil(Mk.iv)) z(ceil(Mk.iv)) z(ceil(Mk.ihp))], ...
                    'color', [1 0 0], 'LineStyle', '-');
                % not significant features
                if Mk.sig == 0
                    for i = 1:length(p)
                        p{i}.FaceAlpha = alpha;
                    end
                    pm.LineStyle = ':';
                end
            end
            % profile
            plot(ax, x, z, 'k', 'LineWidth', linewidth)
            % threshold for "Open" and "Closed"
            if any(Fsig == ["Open", "Closed"]) && ~isnan(NIsig)
                yline(ax, NIsig, '--', 'LineWidth', 1);
            end
        end

        function p = patch_featureelement(~, ax, z, dx, Mr, col)
            dir = sign(Mr.ihp - Mr.ilp);
            ihi = [Mr.ilp; Mr.ihi];
            zlp = z(floor(Mr.ilp));
            p = {};
            for i = 1:2:length(ihi) - 1
                i1 = abs(ceil(dir*ihi(i)));
                i2 = abs(floor(dir*ihi(i+1)));
                xf = dx*([ihi(i); (i1:dir:i2)'; ihi(i+1)] - 1);
                zf = [zlp; z(i1:dir:i2); zlp];
                p{end+1} = patch(ax, [xf; flip(xf)], ...
                    [zf; ones(length(zf), 1)*zlp], col, 'LineStyle', 'none'); %#ok<AGROW>
            end
        end
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % functions of the toolbox (included when compiled)
            if ~isdeployed
                root = fileparts(fileparts(mfilename('fullpath')));
                addpath(fullfile(root, 'featurecharacterization2d'), ...
                    fullfile(root, 'softgauge'));
            end
            app.NIsigEditField.Enable = 'off';
            app.StatisticsvalueEditField.Enable = 'off';
            update_FC(app)
        end

        % Value changed function: FeaturetypeDropDown
        function FeaturetypeDropDownValueChanged(app, event)
            update_FC(app);
            update_plot(app);
        end

        % Button pushed function: LoadprofileButton
        function LoadprofileButtonPushed(app, event)
            [filename, path] = uigetfile({'*.smd;*.mat', ...
                'Profiles (*.smd, *.mat)'});
            if isequal(filename, 0)
                return
            end
            load_profile(app, fullfile(path, filename))
        end

        % Value changed function: PruningTypeDropDown
        function PruningTypeDropDownValueChanged(app, event)
            PT = convertCharsToStrings(app.PruningTypeDropDown.Value);
            if PT == "None"
                app.ThresholdEditField.Enable = 'off';
            else
                app.ThresholdEditField.Enable = 'on';
            end
            update_FC(app);
            update_plot(app)
        end

        % Value changed function: ThresholdEditField
        function ThresholdEditFieldValueChanged(app, event)
            update_FC(app);
            update_plot(app)
        end

        % Value changed function: FeatureCharacterizationEditField
        function FeatureCharacterizationEditFieldValueChanged(app, event)
            FC_complete = convertCharsToStrings(app.FeatureCharacterizationEditField.Value);
            FC = strtrim(split(FC_complete, ";"));
            % pruning
            str = split(strrep(FC(3), "%", " %"));
            str = str(str ~= "");
            app.PruningTypeDropDown.Value = str(1);
            app.ThresholdEditField.Enable = on_off(str(1) ~= "None");
            if length(str) >= 2 && ~isnan(str2double(str(2)))
                app.ThresholdEditField.Value = str2double(str(2));
            end
            app.CheckBox.Value = any(str == "%");
            % significant features
            str = split(strrep(FC(4), "%", " %"));
            str = str(str ~= "");
            app.SignificantfeaturesDropDown.Value = str(1);
            app.NIsigEditField.Enable = on_off(str(1) ~= "All");
            if length(str) >= 2 && ~isnan(str2double(str(2)))
                app.NIsigEditField.Value = str2double(str(2));
            end
            app.CheckBox_2.Value = any(str == "%");
            % attribute type and statistics
            app.AttributetypeDropDown.Value = FC(5);
            str = split(FC(6));
            str = str(str ~= "");
            app.AttributestatisticsDropDown.Value = str(1);
            if length(str) >= 2 && ~isnan(str2double(str(end)))
                app.StatisticsvalueEditField.Value = str2double(str(end));
                app.StatisticsvalueEditField.Enable = 'on';
            else
                app.StatisticsvalueEditField.Enable = 'off';
            end
            update_plot(app)
        end

        % Value changed function: SignificantfeaturesDropDown
        function SignificantfeaturesDropDownValueChanged(app, event)
            Fsig = convertCharsToStrings(app.SignificantfeaturesDropDown.Value);
            if Fsig == "All"
                app.NIsigEditField.Enable = 'off';
            else
                app.NIsigEditField.Enable = 'on';
            end
            update_FC(app)
            update_plot(app)
        end

        % Value changed function: AttributestatisticsDropDown
        function AttributestatisticsDropDownValueChanged(app, event)
            Astats = convertCharsToStrings(app.AttributestatisticsDropDown.Value);
            if Astats == "Perc" || Astats == "Hist"
                app.StatisticsvalueEditField.Enable = 'on';
            else
                app.StatisticsvalueEditField.Enable = 'off';
            end
            update_FC(app)
            update_plot(app)
        end

        % Value changed function: NIsigEditField
        function NIsigEditFieldValueChanged(app, event)
            update_FC(app)
            update_plot(app)
        end

        % Value changed function: CheckBox
        function CheckBoxValueChanged(app, event)
            update_FC(app)
            update_plot(app)
        end

        % Value changed function: CheckBox_2
        function CheckBox_2ValueChanged(app, event)
            update_FC(app)
            update_plot(app)
        end

        % Value changed function: AttributetypeDropDown
        function AttributetypeDropDownValueChanged(app, event)
            update_FC(app)
            update_plot(app)
        end

        % Value changed function: StatisticsvalueEditField
        function StatisticsvalueEditFieldValueChanged(app, event)
            update_FC(app)
            update_plot(app)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1369 492];
            app.UIFigure.Name = 'Feature Characterization';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Feature Characterization')
            xlabel(app.UIAxes, 'profile length / mm')
            ylabel(app.UIAxes, 'profile height / µm')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Box = 'on';
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.Position = [315 8 1035 474];

            % Create WatershedsegmentationPanel
            app.WatershedsegmentationPanel = uipanel(app.UIFigure);
            app.WatershedsegmentationPanel.Title = 'Watershed segmentation';
            app.WatershedsegmentationPanel.Position = [9 293 293 143];

            % Create CheckBox
            app.CheckBox = uicheckbox(app.WatershedsegmentationPanel);
            app.CheckBox.ValueChangedFcn = createCallbackFcn(app, @CheckBoxValueChanged, true);
            app.CheckBox.Text = '%';
            app.CheckBox.Position = [239 13 33 22];
            app.CheckBox.Value = true;

            % Create PruningTypeDropDownLabel
            app.PruningTypeDropDownLabel = uilabel(app.UIFigure);
            app.PruningTypeDropDownLabel.HorizontalAlignment = 'right';
            app.PruningTypeDropDownLabel.Position = [18 342 76 22];
            app.PruningTypeDropDownLabel.Text = 'Pruning Type';

            % Create PruningTypeDropDown
            app.PruningTypeDropDown = uidropdown(app.UIFigure);
            app.PruningTypeDropDown.Items = {'None', 'Wolfprune', 'Width', 'VolS', 'DevLength'};
            app.PruningTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @PruningTypeDropDownValueChanged, true);
            app.PruningTypeDropDown.Position = [128 342 152 22];
            app.PruningTypeDropDown.Value = 'Wolfprune';

            % Create FeaturetypeDropDownLabel
            app.FeaturetypeDropDownLabel = uilabel(app.UIFigure);
            app.FeaturetypeDropDownLabel.HorizontalAlignment = 'right';
            app.FeaturetypeDropDownLabel.Position = [17 378 72 22];
            app.FeaturetypeDropDownLabel.Text = 'Feature type';

            % Create FeaturetypeDropDown
            app.FeaturetypeDropDown = uidropdown(app.UIFigure);
            app.FeaturetypeDropDown.Items = {'H', 'D', 'P', 'V'};
            app.FeaturetypeDropDown.ValueChangedFcn = createCallbackFcn(app, @FeaturetypeDropDownValueChanged, true);
            app.FeaturetypeDropDown.Position = [130 378 150 22];
            app.FeaturetypeDropDown.Value = 'D';

            % Create ThresholdEditFieldLabel
            app.ThresholdEditFieldLabel = uilabel(app.UIFigure);
            app.ThresholdEditFieldLabel.HorizontalAlignment = 'right';
            app.ThresholdEditFieldLabel.Position = [17 306 59 22];
            app.ThresholdEditFieldLabel.Text = 'Threshold';

            % Create ThresholdEditField
            app.ThresholdEditField = uieditfield(app.UIFigure, 'numeric');
            app.ThresholdEditField.ValueChangedFcn = createCallbackFcn(app, @ThresholdEditFieldValueChanged, true);
            app.ThresholdEditField.Position = [131 306 105 22];
            app.ThresholdEditField.Value = 5;

            % Create LoadprofileButton
            app.LoadprofileButton = uibutton(app.UIFigure, 'push');
            app.LoadprofileButton.ButtonPushedFcn = createCallbackFcn(app, @LoadprofileButtonPushed, true);
            app.LoadprofileButton.Position = [9 451 293 29];
            app.LoadprofileButton.Text = 'Load profile';

            % Create FeatureparameterPanel
            app.FeatureparameterPanel = uipanel(app.UIFigure);
            app.FeatureparameterPanel.Title = 'Feature parameter';
            app.FeatureparameterPanel.Position = [11 79 291 207];

            % Create SignificantfeaturesDropDownLabel
            app.SignificantfeaturesDropDownLabel = uilabel(app.FeatureparameterPanel);
            app.SignificantfeaturesDropDownLabel.HorizontalAlignment = 'right';
            app.SignificantfeaturesDropDownLabel.Position = [11 149 107 22];
            app.SignificantfeaturesDropDownLabel.Text = 'Significant features';

            % Create SignificantfeaturesDropDown
            app.SignificantfeaturesDropDown = uidropdown(app.FeatureparameterPanel);
            app.SignificantfeaturesDropDown.Items = {'All', 'Open', 'Closed', 'Top', 'Bot'};
            app.SignificantfeaturesDropDown.ValueChangedFcn = createCallbackFcn(app, @SignificantfeaturesDropDownValueChanged, true);
            app.SignificantfeaturesDropDown.Position = [121 149 148 22];
            app.SignificantfeaturesDropDown.Value = 'All';

            % Create AttributestatisticsDropDownLabel
            app.AttributestatisticsDropDownLabel = uilabel(app.FeatureparameterPanel);
            app.AttributestatisticsDropDownLabel.HorizontalAlignment = 'right';
            app.AttributestatisticsDropDownLabel.Position = [11 44 99 22];
            app.AttributestatisticsDropDownLabel.Text = 'Attribute statistics';

            % Create AttributestatisticsDropDown
            app.AttributestatisticsDropDown = uidropdown(app.FeatureparameterPanel);
            app.AttributestatisticsDropDown.Items = {'Mean', 'Max', 'Min', 'StdDev', 'Perc', 'Hist', 'Sum', 'Density'};
            app.AttributestatisticsDropDown.ValueChangedFcn = createCallbackFcn(app, @AttributestatisticsDropDownValueChanged, true);
            app.AttributestatisticsDropDown.Position = [121 44 148 22];
            app.AttributestatisticsDropDown.Value = 'Mean';

            % Create NIsigEditFieldLabel
            app.NIsigEditFieldLabel = uilabel(app.FeatureparameterPanel);
            app.NIsigEditFieldLabel.HorizontalAlignment = 'right';
            app.NIsigEditFieldLabel.Position = [13 113 32 22];
            app.NIsigEditFieldLabel.Text = 'NIsig';

            % Create NIsigEditField
            app.NIsigEditField = uieditfield(app.FeatureparameterPanel, 'numeric');
            app.NIsigEditField.ValueChangedFcn = createCallbackFcn(app, @NIsigEditFieldValueChanged, true);
            app.NIsigEditField.Position = [122 113 103 22];

            % Create StatisticsvalueEditFieldLabel
            app.StatisticsvalueEditFieldLabel = uilabel(app.FeatureparameterPanel);
            app.StatisticsvalueEditFieldLabel.HorizontalAlignment = 'right';
            app.StatisticsvalueEditFieldLabel.Position = [13 11 85 22];
            app.StatisticsvalueEditFieldLabel.Text = 'Statistics value';

            % Create StatisticsvalueEditField
            app.StatisticsvalueEditField = uieditfield(app.FeatureparameterPanel, 'numeric');
            app.StatisticsvalueEditField.ValueChangedFcn = createCallbackFcn(app, @StatisticsvalueEditFieldValueChanged, true);
            app.StatisticsvalueEditField.Position = [122 11 147 22];

            % Create AttributetypeDropDownLabel
            app.AttributetypeDropDownLabel = uilabel(app.FeatureparameterPanel);
            app.AttributetypeDropDownLabel.HorizontalAlignment = 'right';
            app.AttributetypeDropDownLabel.Position = [12 77 75 22];
            app.AttributetypeDropDownLabel.Text = 'Attribute type';

            % Create AttributetypeDropDown
            app.AttributetypeDropDown = uidropdown(app.FeatureparameterPanel);
            app.AttributetypeDropDown.Items = {'HDh', 'HDw', 'HDv', 'HDl', 'PVh', 'Curvature', 'Count'};
            app.AttributetypeDropDown.ValueChangedFcn = createCallbackFcn(app, @AttributetypeDropDownValueChanged, true);
            app.AttributetypeDropDown.Position = [121 77 149 22];
            app.AttributetypeDropDown.Value = 'HDh';

            % Create CheckBox_2
            app.CheckBox_2 = uicheckbox(app.FeatureparameterPanel);
            app.CheckBox_2.ValueChangedFcn = createCallbackFcn(app, @CheckBox_2ValueChanged, true);
            app.CheckBox_2.Text = '%';
            app.CheckBox_2.Position = [237 113 33 22];

            % Create ConventionISO21920Panel
            app.ConventionISO21920Panel = uipanel(app.UIFigure);
            app.ConventionISO21920Panel.Title = 'Convention ISO-21920';
            app.ConventionISO21920Panel.Position = [12 8 292 60];

            % Create FeatureCharacterizationEditField
            app.FeatureCharacterizationEditField = uieditfield(app.ConventionISO21920Panel, 'text');
            app.FeatureCharacterizationEditField.ValueChangedFcn = createCallbackFcn(app, @FeatureCharacterizationEditFieldValueChanged, true);
            app.FeatureCharacterizationEditField.Position = [11 9 209 22];

            % Create Label
            app.Label = uilabel(app.ConventionISO21920Panel);
            app.Label.FontSize = 13;
            app.Label.FontWeight = 'bold';
            app.Label.Position = [244 9 46 22];
            app.Label.Text = '0';

            % Create Label_2
            app.Label_2 = uilabel(app.ConventionISO21920Panel);
            app.Label_2.Position = [228 9 16 22];
            app.Label_2.Text = '=';

            % Create mLabel
            app.mLabel = uilabel(app.ConventionISO21920Panel);
            app.mLabel.FontSize = 13;
            app.mLabel.FontWeight = 'bold';
            app.mLabel.Position = [267 9 47 22];
            app.mLabel.Text = 'µm';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    methods (Access = public)
        % load a profile (*.smd or *.mat) and update the plot
        function load_profile(app, filepath)
            [~, ~, ext] = fileparts(filepath);
            if strcmpi(ext, '.smd')
                [z_file, ~, ~, dx_file] = smd2mat(filepath);
            else
                profile = load(filepath);
                z_file = profile.z;
                if isfield(profile, 'dx')
                    dx_file = profile.dx;
                else
                    dx_file = 0.5e-3; % default step size in mm
                end
            end
            app.z = z_file(:) - mean(z_file);
            app.dx = dx_file;
            update_plot(app)
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = GUI

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

%% unit of the feature parameter (z in µm, dx in mm)
function u = parameter_unit(AT, Astats)
switch AT
    case {"HDh", "PVh"}
        u = "µm";
    case {"HDw", "HDl"}
        u = "mm";
    case "HDv"
        u = "ml/m²";
    case "Curvature"
        u = "1/µm";
    otherwise
        u = "";
end
switch Astats
    case "Perc"
        u = "";
    case "Density"
        u = u + "/cm";
end
end

%% 'on' or 'off' for the property Enable
function s = on_off(state)
if state
    s = 'on';
else
    s = 'off';
end
end
