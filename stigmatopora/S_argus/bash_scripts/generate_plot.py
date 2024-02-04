#!/usr/bin/env python3
# coding: utf-8
"""

.. module:: generate_plot
   :synopsis: This module produces a graphic summary for BUSCO runs based on short summary files
.. versionadded:: 2.0.0
.. versionchanged:: 4.0.0

 This module produces a graphic summary for BUSCO runs based on short summary files

(``python3 generate_plot.py -h`` and user guide for details on how to do it)

Place the short summary files of all BUSCO runs you would like to see on the figure a single folder.
Keep the file named as follow: ``short_summary.[generic|specific].dataset.label.txt``, 'label' being used in the plot as species name

This tool produces the R code of the figure and uses ggplot2 (2.2.0+). If your system is able to run R, this script
automatically runs it.

You can find both the resulting R script for customisation and the figure in the working directory.

Copyright (c) 2016-2023, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.

"""

import os
import sys
import time
import traceback
import argparse
import subprocess
import glob
from shutil import which
from argparse import RawTextHelpFormatter
import logging
from busco.BuscoLogger import BuscoLogger

#: working directory
_plot_dir = ""
#: r file name
_r_file = "busco_figure.R"

# to avoid running R
_no_r = False

#: Get an instance of _logger for keeping track of events
_logger = BuscoLogger.get_logger(__name__)

RCODE = (
    "######################################\n"
    "#\n"
    "# BUSCO summary figure\n"
    "# @version 4.0.0\n"
    "# @since BUSCO 2.0.0\n"
    "# \n"
    "# Copyright (c) 2016-2023, Evgeny Zdobnov (ez@ezlab.org)\n"
    "# Licensed under the MIT license. See LICENSE.md file.\n"
    "#\n"
    "######################################\n"
    "\n"
    "# Load the required libraries\n"
    "library(ggplot2)\n"
    'library("grid")\n'
    "\n"
    "# !!! CONFIGURE YOUR PLOT HERE !!! \n"
    "# Output\n"
    'my_output <- paste(%s1,"busco_figure.png",sep="/") \n'
    "my_width <- 20\n"
    "my_height <- 15\n"
    'my_unit <- "cm"\n'
    "\n"
