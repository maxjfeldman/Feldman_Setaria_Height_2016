#!/bin/bash

##### sv_z2500 ##### 
# Frame 0
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z2500_1/vis_sv_z2500_L2_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:0,zoom:z2500 \
-T 20 \
-c

PID1=$!
wait $PID1

# Frame 90 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z2500_1/vis_sv_z2500_L2_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:90,zoom:z2500 \
-T 20 

PID2=$!
wait $PID2


# Frame 180 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z2500_1/vis_sv_z2500_L2_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:180,zoom:z2500 \
-T 20 

PID3=$!
wait $PID3


# Frame 270 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z2500_1/vis_sv_z2500_L2_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:270,zoom:z2500 \
-T 20 

PID4=$!
wait $PID4

##### sv_z1000 ##### 

# Frame 0 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z1000_1/vis_sv_z1000_L2_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:0,zoom:z1000 \
-T 20

PID5=$!
wait $PID5

# Frame 90 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z1000_1/vis_sv_z1000_L2_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:90,zoom:z1000 \
-T 20 

PID6=$!
wait $PID6

# Frame 180 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z1000_1/vis_sv_z1000_L2_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:180,zoom:z1000 \
-T 20 

PID7=$!
wait $PID7

# Frame 270 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z1000_1/vis_sv_z1000_L2_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:270,zoom:z1000 \
-T 20 

PID8=$!
wait $PID8


##### sv_z500 ##### 

# Frame 0 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z500_1/vis_sv_z500_L1_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:0,zoom:z500 \
-T 20 

PID9=$!
wait $PID9

# Frame 90 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z500_1/vis_sv_z500_L1_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:90,zoom:z500 \
-T 20 

PID10=$!
wait $PID10


# Frame 180 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z500_1/vis_sv_z500_L1_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:180,zoom:z500 \
-T 20 

PID11=$!
wait $PID11


# Frame 270 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z500_1/vis_sv_z500_L1_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:270,zoom:z500 \
-T 20 

PID12=$!
wait $PID12


##### sv_z1_h1 ##### 

# Frame 0 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z1_1/vis_sv_z1_L1_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:0,zoom:z1 \
-D 2014-02-08-05-00-00-00_2014-02-11-03-00-00-00 \
-T 20 

PID13=$!
wait $PID13

# Frame 90 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z1_1/vis_sv_z1_L1_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:90,zoom:z1 \
-D 2014-02-08-05-00-00-00_2014-02-11-03-00-00-00 \
-T 20 

PID14=$!
wait $PID14


# Frame 180 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z1_1/vis_sv_z1_L1_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:180,zoom:z1 \
-D 2014-02-08-05-00-00-00_2014-02-11-03-00-00-00 \
-T 20 

PID15=$!
wait $PID15


# Frame 270 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z1_1/vis_sv_z1_L1_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:270,zoom:z1 \
-D 2014-02-08-05-00-00-00_2014-02-11-03-00-00-00 \
-T 20 

PID16=$!
wait $PID16


##### sv_z1_h2 ##### 

# Frame 0 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z1_2/vis_sv_z1_L0_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:0,zoom:z1 \
-D 2014-02-11-05-00-00-00_2014-02-16-00-00-00-00 \
-T 20 \

PID17=$!
wait $PID17

# Frame 90 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z1_2/vis_sv_z1_L0_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:90,zoom:z1 \
-D 2014-02-11-05-00-00-00_2014-02-16-00-00-00-00 \
-T 20 

PID18=$!
wait $PID18


# Frame 180 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z1_2/vis_sv_z1_L0_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:180,zoom:z1 \
-D 2014-02-11-05-00-00-00_2014-02-16-00-00-00-00 \
-T 20 

PID19=$!
wait $PID19


# Frame 270 
time /home/mfeldman/ril_test_img/ril_out/image_analysis_MF.py \
-d /home/nfahlgren/projects/lemnatec/setaria_exp1/images \
-p /home/mfeldman/ril_test_img/sv_z1_2/vis_sv_z1_L0_e82.py \
-s /home/mfeldman/ril_test_img/ril_out/setaria_ril_db.sqlite3 \
-i /home/mfeldman/ril_test_img/ril_out/ril_output \
-f imgtype_camera_frame_zoom_other \
-M imgtype:VIS,camera:SV,frame:270,zoom:z1 \
-D 2014-02-11-05-00-00-00_2014-02-16-00-00-00-00 \
-T 20 

PID20=$!
wait $PID20