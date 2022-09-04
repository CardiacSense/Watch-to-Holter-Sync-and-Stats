[file, folder] = uigetfile('*.XML'); %Choose an XML file
% Handles some directory structure pecularities in some systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fslash = strfind(folder,'/');
if isempty(fslash)==1
    dirdemarc = '\';
else
    dirdemarc = '/';
end

slashcheck = strcmp(folder(end),dirdemarc);
if slashcheck == 0
    folder = [folder dirdemarc];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%CAW: Attempts to use xmlread failed under Octave, which appears
% to be a well known issue in the community. For MATLAB users
% this script can be made much more robust by using xmlread()
% and the DOM functions to search for the <parsedwaveforms>
% element.

%{
Jack - Try something like this:
xdoc = xmlread('test.xml');
header = xdoc.getElementsByTagName('Header');
id = char(header.item(0).getAttribute('id'));
com = xdoc.getElementsByTagName('Communication');
sub = com.item(0).getElementsByTagName('SubNetwork');
subname = char(sub.item(0).getAttribute('name'));
produces:
id =
IED15
subname =
Subnet2
%}

% CAW: Below you'll find both Octave and Matlab code to read XML
% Comment out the one you don't need, and uncomment the one you do need.

%%
%% Octave XML
%%
% javaaddpath ('D:/Users/Christopher/Downloads/xerces-2_12_1/xercesImpl.jar');
% javaaddpath ('D:/Users/Christopher/Downloads/xerces-2_12_1/xml-apis.jar');
% 
% parser = javaObject("org.apache.xerces.parsers.DOMParser");
% parser.parse([folder file]); 
% xdoc = parser.getDocument();
%%
%% END Octave XML
%%

%%
%% Matlab XML
%%

xdoc = xmlread([folder file]);

%%
%% END Matlab XML
%%

channels = xdoc.getElementsByTagName('CHANNEL');

% Extract each of the 12 leads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for n = 1:12

    % 0-indexing for XML items
    channelXml = channels.item(n-1);
    
    % Divide each sample by this value to get mV
    unitsPerMv = str2num(channelXml.getAttribute('UNITS_PER_MV'));
    
    % Sampling frequency (Hz)
    freq = str2num(channelXml.getAttribute('SAMPLE_FREQ'));
    
    % Number of samples
    samples = str2num(channelXml.getAttribute('DURATION'));
    
    % Decode the base64 encoded data
    b64data = regexprep(char(channelXml.getAttribute('DATA')), '\s+', '');
    data = uint8(base64decode(b64data));
    % But only keep 2*samples bytes worth (sometimes there is a trailing 0)
    data = data(1:(2*samples));
    
    % The samples are 16-bits, so reinterpret them as int16 then convert to doubles
    output = double(typecast(data, 'int16'));

    % Add the next lead to the Cell array
    leads(n) = output / unitsPerMv;
end

%%%
% CAW: I could not get this to work on Octave...
% So I just did the cell2mat calls below in the plotting routines.
%
% Convert leads cell array to numeric matrix
%leads = cell2mat(leads);
%%%

% 12-Lead Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('numbertitle','off','name','12-Lead Plots');
subplot(6,2,1);
plot(cell2mat(leads(1)));
title('I');

subplot(6,2,2);
plot(cell2mat(leads(2)));
title('II');

subplot(6,2,3);
plot(cell2mat(leads(3)));
title('III');

subplot(6,2,4);
plot(cell2mat(leads(4)));
title('aVR');

subplot(6,2,5);
plot(cell2mat(leads(5)));
title('aVL');

subplot(6,2,6);
plot(cell2mat(leads(6)));
title('aVF');

subplot(6,2,7);
plot(cell2mat(leads(7)));
title('V1');

subplot(6,2,8);
plot(cell2mat(leads(8)));
title('V2');

subplot(6,2,9);
plot(cell2mat(leads(9)));
title('V3');

subplot(6,2,10);
plot(cell2mat(leads(10)));
title('V4');

subplot(6,2,11);
plot(cell2mat(leads(11)));
title('V5');

subplot(6,2,12);
plot(cell2mat(leads(12)));
title('V6');