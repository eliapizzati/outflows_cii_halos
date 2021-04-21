# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 15:21:19 2021

@author: anna
"""

import os
import numpy as np
#import matplotlib.pyplot as plt

import natconst as nc


script_dir = os.getcwd()

main_dir = os.path.dirname(script_dir)

data_dir = os.path.join(main_dir,"data")

log_dir = os.path.join(main_dir, "log")

plot_dir = os.path.join(main_dir,"plot")

if not os.path.exists(data_dir):
    os.mkdir(data_dir)

if not os.path.exists(log_dir):
    os.mkdir(log_dir)

if not os.path.exists(plot_dir):
    os.mkdir(plot_dir)


        
