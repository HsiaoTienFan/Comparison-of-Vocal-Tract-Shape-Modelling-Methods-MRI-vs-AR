function data = ext_data(Data, type, i)

if type == 0
data = Data(1).data(: ,i);
end

if type == 1
if i < 3
    if i == 1
  data =  Data(:,i)/10;
    else
  %    data =  Data(:,i)/100;     
 data =  sqrt(Data(:,i)/pi);
    end
else
    data = 0;
end
end