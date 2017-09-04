function tubes = get_tubes_area(data_set_x, data_set_y, tubelength1, start, finish, type)




%set up tube data set 1
if type == 0

    

tube_data = (data_set_x(start)):tubelength1/10:(data_set_x(finish));
tube_1_points = interp1(data_set_x,data_set_y,tube_data);

tube1 = trapz(tube_data(1:11),tube_1_points(1:11))/tubelength1;
tube2 = trapz(tube_data(11:21),tube_1_points(11:21))/tubelength1;
tube3 = trapz(tube_data(21:31),tube_1_points(21:31))/tubelength1;
tube4 = trapz(tube_data(31:41),tube_1_points(31:41))/tubelength1;
tube5 = trapz(tube_data(41:51),tube_1_points(41:51))/tubelength1;
tube6 = trapz(tube_data(51:61),tube_1_points(51:61))/tubelength1;
tube7 = trapz(tube_data(61:71),tube_1_points(61:71))/tubelength1;
tube8 = trapz(tube_data(71:81),tube_1_points(71:81))/tubelength1;
tube9 = trapz(tube_data(81:91),tube_1_points(81:91))/tubelength1;
tube10 = trapz(tube_data(91:101) ,tube_1_points(91:101))/tubelength1;
tube11 = trapz(tube_data(101:111),tube_1_points(101:111))/tubelength1;

else   
    
    
tube_data = (data_set_x(1)):tubelength1/10:(data_set_x(size(data_set_x),1));
tube_1_points = interp1(data_set_x,data_set_y,tube_data);
tube1 = trapz(tube_data(1:11),tube_1_points(1:11))/tubelength1;
tube2 = trapz(tube_data(11:21),tube_1_points(11:21))/tubelength1;
tube3 = trapz(tube_data(21:31),tube_1_points(21:31))/tubelength1;
tube4 = trapz(tube_data(31:41),tube_1_points(31:41))/tubelength1;
tube5 = trapz(tube_data(41:51),tube_1_points(41:51))/tubelength1;
tube6 = trapz(tube_data(51:61),tube_1_points(51:61))/tubelength1;
tube7 = trapz(tube_data(61:71),tube_1_points(61:71))/tubelength1;
tube8 = trapz(tube_data(71:81),tube_1_points(71:81))/tubelength1;
tube9 = trapz(tube_data(81:91),tube_1_points(81:91))/tubelength1;
tube10 = trapz(tube_data(91:101) ,tube_1_points(91:101))/tubelength1;
tube11 = trapz(tube_data(101:111),tube_1_points(101:111))/tubelength1;

end
tubes = [sqrt(tube1/pi), sqrt(tube2/pi), sqrt(tube3/pi), sqrt(tube4/pi), sqrt(tube5/pi), sqrt(tube6/pi), sqrt(tube7/pi), sqrt(tube8/pi), sqrt(tube9/pi), sqrt(tube10/pi), sqrt(tube11/pi), 0];
