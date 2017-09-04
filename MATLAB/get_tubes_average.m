function tubes = get_tubes_average(data_set_x, data_set_y, tubelength1, start, finish, type)




%set up tube data set 1
if type == 0

    

tube_data = (data_set_x(start)):tubelength1/10:(data_set_x(finish));
tube_1_points = interp1(data_set_x,data_set_y,tube_data);

tube1 = mean(tube_1_points(1:11));
tube2 = mean(tube_1_points(11:21));
tube3 = mean(tube_1_points(21:31));
tube4 = mean(tube_1_points(31:41));
tube5 = mean(tube_1_points(41:51));
tube6 = mean(tube_1_points(51:61));
tube7 = mean(tube_1_points(61:71));
tube8 = mean(tube_1_points(71:81));
tube9 = mean(tube_1_points(81:91));
tube10 = mean(tube_1_points(91:101));
tube11 = mean(tube_1_points(101:111));

else   
    
    
tube_data = (data_set_x(1)):tubelength1/10:(data_set_x(size(data_set_x),1));
tube_1_points = interp1(data_set_x,data_set_y,tube_data);

tube1 = mean(tube_1_points(1:11))/2;
tube2 = mean(tube_1_points(11:21))/2;
tube3 = mean(tube_1_points(21:31))/2;
tube4 = mean(tube_1_points(31:41))/2;
tube5 = mean(tube_1_points(41:51))/2;
tube6 = mean(tube_1_points(51:61))/2;
tube7 = mean(tube_1_points(61:71))/2;
tube8 = mean(tube_1_points(71:81))/2;
tube9 = mean(tube_1_points(81:91))/2;
tube10 = mean(tube_1_points(91:101))/2;
tube11 = mean(tube_1_points(101:111))/2;

end
tubes = [sqrt(tube1/pi), sqrt(tube2/pi), sqrt(tube3/pi), sqrt(tube4/pi), sqrt(tube5/pi), sqrt(tube6/pi), sqrt(tube7/pi), sqrt(tube8/pi), sqrt(tube9/pi), sqrt(tube10/pi), sqrt(tube11/pi)];
