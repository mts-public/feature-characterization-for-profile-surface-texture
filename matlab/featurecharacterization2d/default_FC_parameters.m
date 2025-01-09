function xFC = default_FC_parameters(z, dx)
% This function calculates all named parameters of Feature Characterization 
% defined by ISO 21920-2 with the default settings according to ISO 21920-3

%% Watershed segmentation with default pruning setting
Mp = watershed_segmentation(z, dx, "P", "Wolfprune", (5/100)*Rz(z, dx));
Mv = watershed_segmentation(z, dx, "V", "Wolfprune", (5/100)*Rz(z, dx));
%% default feature parameters according to ISO 21920-2
xFC.Rpd = feature_parameter(z, dx, Mp, "All", 1 ,"Count", "Density");
xFC.Rvd = feature_parameter(z, dx, Mv, "All", 1 ,"Count", "Density");
xFC.Rmpc = feature_parameter(z, dx, Mp, "All", 1 ,"Curvature", "Mean");
xFC.Rmvc = feature_parameter(z, dx, Mv, "All", 1 ,"Curvature", "Mean");
xFC.R5p = feature_parameter(z, dx, Mp, "Top", 5,"PVh", "Mean");
xFC.R5v = feature_parameter(z, dx, Mv, "Bot", 5,"PVh", "Mean");
xFC.R10z = xFC.R5p + xFC.R5v;
end
