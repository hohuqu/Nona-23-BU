import streamlit as st
import pandas as pd
import numpy as np
import os
from Bio import SeqIO
from io import StringIO

# Note: to run command 'streamlit run *.py', make sure you are at the correct
# directory where the file is

st.title('BioGenius') # Name of the App

uploaded_file = st.file_uploader("Upload your file here...", type=['fasta'])

if uploaded_file is not None:
    
    fasta_content = uploaded_file.read().decode("utf-8")
    st.write(f"{uploaded_file.name}:")
    
    sequences = list(SeqIO.parse(StringIO(fasta_content), "fasta"))
    for record in sequences:
        st.write(f">{record.id}")
        st.write(record.seq)
        
    # PCA:
    st.markdown('## <span style="font-size:30px;">Polymerase Cycling Assembly (PCA)</span>', unsafe_allow_html=True)
    for record in sequences:
        st.write(f">{record.id}")
        # output = record.myFunction(record.seq)
        output1 = "List generated genes here" # switch to function
        st.write(output1)
        output2 = "List the PCA reactions that would work" # switch to function
        st.write(output2)
        output3 = "List the possible recombinant DNA molecules" # switch to function
        st.write(output3)
      
    # Restriction
    st.markdown('## <span style="font-size:30px;">Restriction Assembly</span>', unsafe_allow_html=True)
    for record in sequences:
        st.write(f">{record.id}")
        # output = record.myFunction(record.seq)
        output1 = "Ouput 1" # switch to function
        st.write(output1)
        output2 = "Ouput 2" # switch to function
        st.write(output2)
        output3 = "Ouput 3" # switch to function
        st.write(output3)
    
    # Gibson Assesmbly
    st.markdown('## <span style="font-size:30px;">Gibson Assembly</span>', unsafe_allow_html=True)
    for record in sequences:
        st.write(f">{record.id}")
        # output = record.myFunction(record.seq)
        output1 = "Ouput 1" # switch to function
        st.write(output1)
        output2 = "Ouput 2" # switch to function
        st.write(output2)
        output3 = "Ouput 3" # switch to function
        st.write(output3)
    
