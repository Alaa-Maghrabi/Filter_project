### This folder contains everything you need to know about experiment A ###

You run on matlab the file Experiment.m
and you load any of the 24 experiments located on the same folder.
The experiments called humanX are with human around.
The experiments called humanless are with no human around (I switch the tag's 
position between the side of my body and in front of my body)


This file is composed of two parts:

first part:
   Run three filters: kalman filter, Standard Particle filter, Particle filter using
kalman's velocity
-The kalman filter is implemented in the file kalman.m
-The two other filters are implemented in the file particle_filter.m: you need to 
uncomment the different models to switch between the standard particle filter and 
the one using kalman's velocity


Second part:
  Contains only the Kalman-Like particle filter called particle_fil.m which calls a
kalman filter named kalm.m
