#!/bin/bash	

multiplet=2 
user_name=smph-kitkat # need to be found according to user id
user_id=1510 # to be retrieved from timestamp file
currentDate=$(date +%s)

g++ analysis_main.cpp `root-config --cflags` `root-config --glibs` -o analysisMain

./analysisMain $multiplet $user_name $user_id $currentDate
#./analysisMain $multiplet $user_name $user_id $currentDate>output.txt  ---> uncomment this line if you want to
																			#pipe the output in a txt file  
root -l "plot4user.C(\"$user_name\")"