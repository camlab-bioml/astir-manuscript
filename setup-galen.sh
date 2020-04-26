#!/bin/bash

# This bash script creates data-raw & corresponding subdirectories in the current directory
# and creates symbolic links from mass cytometry data stored in /home/ltri/campbell/share/datasets/
# to be used as raw data input

mkdir -p data-raw
mkdir -p data-raw/wagner-2019

ln -s /home/ltri/campbell/share/datasets/jackson-imc/ data-raw/
ln -s /home/ltri/campbell/share/datasets/wagner-2019-cell/fcs_files data-raw/wagner-2019/fcs_files

