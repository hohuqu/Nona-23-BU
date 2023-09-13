# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 14:54:19 2023

@author: Yehuda and Jade
"""

import streamlit as st
import pandas as pd
import numpy as np
import os
from Bio import SeqIO
from io import StringIO

# Note: to run command 'streamlit run *.py', make sure you are at the correct
# directory where the file is

st.title('Fasta Input')

uploaded_file = st.file_uploader("Upload your file here...", type=['fasta'])

if uploaded_file is not None:
    
    fasta_content = uploaded_file.read().decode("utf-8")
    st.write(f"Reading and parsing {uploaded_file.name}...")
    
    sequences = SeqIO.parse(StringIO(fasta_content), "fasta")
    for record in sequences:
        st.write(f">{record.id}")
        st.write(record.seq)
    
