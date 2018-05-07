#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 10:06:10 2018

@author: bioinformatics
"""


# WRITE TO FILE
import pickle 
file_results = open('gnomex_paired_end_single_index.obj', 'w') 
pickle.dump(results, file_results) 
close(file_results)

# file_results = open('gnomex_paired_end_dual_index.obj', 'w') 
# pickle.dump(results, file_results) 
# close(file_results)

# IMPORT TO PYTHON 
 fh = open('gnomex_paired_end_dual_index.obj', 'r')
 results = pickle.load(fh)