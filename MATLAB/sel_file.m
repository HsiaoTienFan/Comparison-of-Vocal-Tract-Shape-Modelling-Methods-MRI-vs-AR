function Data = sel_file(datatype)

if datatype == 0
[filename, pathname] = uigetfile({'*.000', 'R-data (*.000)';'*.*',  'All Files (*.*)'}, 'Select a AR data file');
current = cd;
cd(pathname);
files = dir('*.000');
for i = 1:length(files)
Data(i) = importdata(files(i).name);

end
cd(current)
end
if datatype == 1
        
[filename, pathname] = uigetfile({'*.txt', 'Text (*.txt)';'*.*',  'All Files (*.*)'}, 'Select a MRI data file');
current = cd;
cd(pathname);
files = dir('*.txt');
for i = 1:length(files)
Data{i} = importdata(files(i).name);

end
cd(current)









 
end



    