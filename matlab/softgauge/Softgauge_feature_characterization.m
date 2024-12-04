clc;clear

[file,path] = uigetfile("*.smd");
filepath = append(path, file);
[z, ~, ~, dx] = smd2mat(filepath);

xFC = default_FC_parameters(z, dx);