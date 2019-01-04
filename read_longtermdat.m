function [C_data, C_struct, C_header] = read_longtermdat

% "day" "blk" "trl" "stim" "first" "second" "RT" "nstim" "name" "delta" "laststim" "jump" "lastjump" "cor1" "cor2"
% "1" "blk00" "trl000" 13 9 12 3587 30 "Monique" -4 NA NA NA 0 0
% NOTE: Last two columns were originally TRUE, FALSE. Changed to 0, 1.

debug_mode = 'off';

switch debug_mode
    case 'off'
        filename = 'fatih_longterm.dat';
    case 'on'
        filename = 'toy_longterm.dat';
end

fid = fopen(filename);
formatSpec = '%q';
N = 15;
C_header = textscan(fid, formatSpec, N, 'Delimiter',' ');
C_data = textscan(fid, ['%q %q %q %n %n %n %n %n %q %n %n ' ...
    '%n %n %n %n'], 'TreatAsEmpty', {'NA','na'});
fclose(fid);

% Data size
n = length(C_data{1});

% Put all in a structure
for i = n:-1:1
    C_struct(i).day = str2double(C_data{1}{i});
    C_struct(i).blk = C_data{2}{i};
    C_struct(i).trl = C_data{3}{i};
    C_struct(i).stim = C_data{4}(i);
    C_struct(i).first = C_data{5}(i);
    C_struct(i).second = C_data{6}(i);
    C_struct(i).RT = C_data{7}(i);
    C_struct(i).nstim = C_data{8}(i);
    C_struct(i).name = C_data{9}{i};
    C_struct(i).delta = C_data{10}(i);
    C_struct(i).lastsim = C_data{11}(i);
    C_struct(i).jump = C_data{12}(i);
    C_struct(i).lastjump = C_data{13}(i);
    C_struct(i).cor1 = C_data{14}(i);
    C_struct(i).cor2 = C_data{15}(i);
end
    
end