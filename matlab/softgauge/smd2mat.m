function [z, L, x, dx] = smd2mat(filepath)
% reads a profile from a softgauge file (*.smd) according to ISO 5436-2.
% The file consists of four sections separated by ETX (char(3)):
%   1. header with the axis definitions, e.g.
%        CX I 8001 mm 1.0E+000 D 5.0E-004  (incremental x-axis, increment)
%        CZ A 8001 um 1.0E+000 D           (absolute z-axis)
%      fields: axis, axis type, number of points, unit, scale factor,
%      data type (and increment for incremental axes)
%   2. additional information (e.g. DATE, CREATED-BY)
%   3. z-values, one value per line
%   4. checksum: sum of all bytes up to the end of the line of the third
%      ETX modulo 65535
% Fields can be separated by NUL characters or blanks.
% INPUTS:
%   filepath - path of the *.smd file
% OUTPUTS:
%   z        - vertical profile values in µm (column vector)
%   L        - profile length n*dx in mm
%   x        - x-positions in mm (row vector)
%   dx       - step size in x-direction in mm

fid = fopen(filepath, 'r');
if fid < 0
    error("File '%s' cannot be opened.", filepath)
end
bytes = fread(fid, Inf, '*uint8')';
fclose(fid);

%% sections
iETX = find(bytes == 3);
if length(iETX) < 3
    error("'%s' is not a valid smd file (less than 3 ETX separators).", ...
        filepath)
end
header = char(bytes(1:iETX(1)-1));
data = char(bytes(iETX(2)+1:iETX(3)-1));

%% header: axis definitions
header(header == char(0)) = ' ';
lines = splitlines(string(header));
CX = []; CZ = [];
for i = 1:length(lines)
    fields = split(strtrim(lines(i)));
    switch fields(1)
        case "CX"
            CX = fields;
        case "CZ"
            CZ = fields;
    end
end
if isempty(CX) || isempty(CZ)
    error("'%s' does not contain a CX and a CZ axis definition.", filepath)
end
if CX(2) ~= "I" || length(CX) < 7
    error("Only incremental x-axes (CX I ... increment) are supported.")
end
if CZ(2) ~= "A"
    error("Only absolute z-axes (CZ A ...) are supported.")
end
dx = str2double(CX(7))*str2double(CX(5))*unit_factor(CX(4), "mm");
scale_z = str2double(CZ(5))*unit_factor(CZ(4), "um");

%% z-values
lines = splitlines(strtrim(string(data)));
z = str2double(strtrim(lines));
if any(isnan(z))
    error("Invalid z-value in line %d of the data section.", ...
        find(isnan(z), 1))
end
z = z*scale_z;
n = length(z);
if n ~= str2double(CZ(3))
    warning("Number of z-values (%d) differs from the header (%s).", ...
        n, CZ(3))
end

%% checksum
if length(iETX) >= 4
    % end of the line of the third ETX
    iEOL = iETX(3) + find(bytes(iETX(3)+1:iETX(4)) == 10, 1);
    if isempty(iEOL)
        iEOL = iETX(3);
    end
    checksum = str2double(strtrim(string(char(bytes(iEOL+1:iETX(4)-1)))));
    if mod(sum(double(bytes(1:iEOL))), 65535) ~= checksum
        warning("Checksum of '%s' does not match.", filepath)
    end
end

%% x-values
L = n*dx;
x = (0:n-1)*dx;
end

%% conversion factor of a length unit to the target unit
function f = unit_factor(unit, target)
units = ["m", "mm", "um", "nm"];
exponents = [0, -3, -6, -9];
i = find(units == unit, 1);
if isempty(i)
    error("Unknown unit '%s'.", unit)
end
f = 10^(exponents(i) - exponents(units == target));
end
