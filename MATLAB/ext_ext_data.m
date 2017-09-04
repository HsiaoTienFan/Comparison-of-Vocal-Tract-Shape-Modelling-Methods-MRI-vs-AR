function data = ext_ext_data(Data, start, finish, type)

if type == 0
data =Data(start:finish ,1);
end

if type == 1
if i < 3
  data =  Data;
else
    data = 0;
end
end