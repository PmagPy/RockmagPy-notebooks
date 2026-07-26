04/15/2026
FORCme version 1

# FORCme - First Order Reversal Curve, machine-code-execution

# Requires: FORCme.ipynb and forc_funcs.py

# ========================= #
# Project information
# ========================= #

# Code implementation: Maxwell Brown (Institute for Rock Magnetism, University of Minnesota)
# Version: V1 (beta)
# Date: 04/15/2026
# Based on ideas, mathematics, and code developed by Chris Pike, Andrew Roberts, Richard Harrison, Ramon Egli, Adrian Muxworthy (and others). 
# Part of the RockmagPy coding suite (https://pmagpy.github.io/RockmagPy-notebooks/book/intro.html).
# OpenAI ChatGPT 5.2 assisted in executing the mechanics of the coding.
# 
# Work carried out under NSF Award EAR-2148549.
#
# Please cite as:
# Maxwell Brown (2026), FORCme: Application of machine code and adaptive smoothing within python for the rapid processing of First Order Reversal Curves, xxx

# =========
# Setup
# =========

Package contains FORCme Jupyter notebook and forc_funcs.py file.
Add forc_funcs.py file to your PmagPy/pmagpy folder (will also run with forc_funcs.py in the same directory as the jupyter notebook.
Make sure PmagPy/pmagpy is on your python PATH.
Jupyter notebook can be saved anywhere and will call forc_funcs.py if the PATH is correctly setup.

Before running install NUMBA
