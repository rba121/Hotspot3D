open_project hotspot3D.prj -reset
set_top computeTempFPGA

add_files hotspot3D.cpp
add_files -tb hotspot3D_test.cpp

open_solution "solution1"
set_part {xcu50-fsvh2104-2-e}
create_clock -period 3.34

csim_design
csynth_design

close_project
exit
