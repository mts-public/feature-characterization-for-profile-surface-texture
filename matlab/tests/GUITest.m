classdef GUITest < matlab.uitest.TestCase
% Tests of the GUI (matlab/GUI/GUI.m): the displayed feature parameter has
% to match the result of feature_characterization of the toolbox.
% Run from any folder:
%   matlab -batch "results = runtests('matlab/tests/GUITest.m')"

    properties
        App
        Root
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            testCase.Root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(fullfile(testCase.Root, 'matlab', 'GUI'), ...
                fullfile(testCase.Root, 'matlab', 'featurecharacterization2d'), ...
                fullfile(testCase.Root, 'matlab', 'softgauge'));
        end
    end

    methods (TestMethodSetup)
        function launchApp(testCase)
            testCase.App = GUI;
            testCase.addTeardown(@delete, testCase.App);
            testCase.App.load_profile(fullfile(testCase.Root, 'data', ...
                'profiles', 'Bu_1_56_ak.smd'));
        end
    end

    properties (TestParameter)
        config = {
            "FC; D; Wolfprune 5 %; All; HDh; Mean"
            "FC; P; Wolfprune 5 %; Top 5; PVh; Mean"
            "FC; V; Wolfprune 5 %; Bot 5; PVh; Mean"
            "FC; D; Width 5 %; All; HDw; Mean"
            "FC; D; VolS 1; All; HDv; Sum"
            "FC; D; DevLength 0.05; All; HDl; Max"
            "FC; H; Wolfprune 5 %; Closed 50 %; HDh; Mean"
            "FC; D; Wolfprune 5 %; Open 0; HDh; Mean"
            "FC; P; Wolfprune 5 %; All; Curvature; Mean"
            "FC; P; Wolfprune 5 %; All; Count; Density"
            "FC; D; Wolfprune 5 %; All; HDh; Perc 2"
            "FC; D; Wolfprune opt; All; HDh; Mean"
            }
    end

    methods (Test)
        function displayedValueMatchesToolbox(testCase, config)
            testCase.type(testCase.App.FeatureCharacterizationEditField, config);
            [z, ~, ~, dx] = smd2mat(fullfile(testCase.Root, 'data', ...
                'profiles', 'Bu_1_56_ak.smd'));
            FC = strtrim(split(config, ";"));
            xFC = feature_characterization(z - mean(z), dx, FC(2), FC(3), ...
                FC(4), FC(5), FC(6));
            testCase.verifyEqual(string(testCase.App.Label.Text), ...
                string(num2str(xFC, 2)));
            % profile and motifs are plotted
            testCase.verifyGreaterThan(numel(testCase.App.UIAxes.Children), 1);
        end

        function dropDownsUpdateConvention(testCase)
            testCase.choose(testCase.App.FeaturetypeDropDown, 'P');
            testCase.choose(testCase.App.SignificantfeaturesDropDown, 'Top');
            testCase.type(testCase.App.NIsigEditField, 5);
            testCase.choose(testCase.App.AttributetypeDropDown, 'PVh');
            testCase.verifyEqual(string(testCase.App.FeatureCharacterizationEditField.Value), ...
                "FC; P; Wolfprune 5 %; Top 5; PVh; Mean");
            testCase.verifyEqual(string(testCase.App.mLabel.Text), "µm");
        end

        function thresholdLineForOpen(testCase)
            testCase.type(testCase.App.FeatureCharacterizationEditField, ...
                "FC; D; Wolfprune 5 %; Open 0.5; HDh; Mean");
            lines = findobj(testCase.App.UIAxes, 'Type', 'ConstantLine');
            testCase.verifyNumElements(lines, 1);
            testCase.verifyEqual(lines.Value, 0.5);
        end

        function noMotifs(testCase)
            testCase.type(testCase.App.FeatureCharacterizationEditField, ...
                "FC; D; Wolfprune 1000; All; HDh; Mean");
            testCase.verifyEqual(string(testCase.App.Label.Text), "NaN");
        end

        function loadMatFile(testCase)
            testCase.App.load_profile(fullfile(testCase.Root, 'data', ...
                'profiles', 'Bu_1_56_ak.mat'));
            S = load(fullfile(testCase.Root, 'data', 'profiles', 'Bu_1_56_ak.mat'));
            xFC = feature_characterization(S.z - mean(S.z), 0.5e-3, "D", ...
                "Wolfprune 5 %", "All", "HDh", "Mean");
            testCase.verifyEqual(string(testCase.App.Label.Text), ...
                string(num2str(xFC, 2)));
        end
    end
end
