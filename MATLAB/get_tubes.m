function tubes = get_tubes(data_set_x, data_set_y, tubelength1, start, finish, type)

%set up tube data set 1
if type == 0
tube1 = interp1(data_set_x,data_set_y,((tubelength1/2)+data_set_x(start)));
tube2 = interp1(data_set_x,data_set_y,((3*tubelength1/2)+data_set_x(start)));
tube3 = interp1(data_set_x,data_set_y,((5*tubelength1/2)+data_set_x(start)));
tube4 = interp1(data_set_x,data_set_y,((7*tubelength1/2)+data_set_x(start)));
tube5 = interp1(data_set_x,data_set_y,((9*tubelength1/2)+data_set_x(start)));
tube6 = interp1(data_set_x,data_set_y,((11*tubelength1/2)+data_set_x(start)));
tube7 = interp1(data_set_x,data_set_y,((13*tubelength1/2)+data_set_x(start)));
tube8 = interp1(data_set_x,data_set_y,((15*tubelength1/2)+data_set_x(start)));
tube9 = interp1(data_set_x,data_set_y,((17*tubelength1/2)+data_set_x(start)));
tube10 = interp1(data_set_x,data_set_y,((19*tubelength1/2)+data_set_x(start)));
tube11 = interp1(data_set_x,data_set_y,((21*tubelength1/2)+data_set_x(start)));

else 
tube1 = interp1(data_set_x,data_set_y,((tubelength1/2)));
tube2 = interp1(data_set_x,data_set_y,((3*tubelength1/2)));
tube3 = interp1(data_set_x,data_set_y,((5*tubelength1/2)));
tube4 = interp1(data_set_x,data_set_y,((7*tubelength1/2)));
tube5 = interp1(data_set_x,data_set_y,((9*tubelength1/2)));
tube6 = interp1(data_set_x,data_set_y,((11*tubelength1/2)));
tube7 = interp1(data_set_x,data_set_y,((13*tubelength1/2)));
tube8 = interp1(data_set_x,data_set_y,((15*tubelength1/2)));
tube9 = interp1(data_set_x,data_set_y,((17*tubelength1/2)));
tube10 = interp1(data_set_x,data_set_y,((19*tubelength1/2)));
tube11 = interp1(data_set_x,data_set_y,((21*tubelength1/2)));

end
tubes = [sqrt(tube1/pi), sqrt(tube2/pi), sqrt(tube3/pi), sqrt(tube4/pi), sqrt(tube5/pi), sqrt(tube6/pi), sqrt(tube7/pi), sqrt(tube8/pi), sqrt(tube9/pi), sqrt(tube10/pi), sqrt(tube11/pi), 0];
