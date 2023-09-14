import streamlit as st
import pandas as pd
import numpy as np
import os
import pca_gibson
import restriction_site
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO

# Note: to run command 'streamlit run *.py', make sure you are at the correct directory where the file is

# Possible Target Sequences
st.title('BioGenius') # Name of the App

# Choose an assembly
option = st.selectbox(
    'Choose an Assembly:',
    ('Polymerase Cycling Assembly (PCA)', 'Restriction Assembly', 'Gibson Assembly'))

# Upload or drag file
uploaded_file = st.file_uploader("Upload your file here:", type=['fasta'])

if uploaded_file is not None:
    
    fasta_content = uploaded_file.read().decode("utf-8")
    st.write(f"{uploaded_file.name}:")
    
    # Read the contents of the
    sequences = list(SeqIO.parse(StringIO(fasta_content), "fasta"))
    for record in sequences:
        st.write(f">{record.id}")
        st.write(record.seq)
    
    # Create option boxes for selecting 2 oligo sequences
    oligos_options = st.multiselect('Choose Exact 2 Sequences:',[record.id for record in sequences])
    
    if len(oligos_options ) == 2:
        # Display 2 selected sequences
        selected_sequences = [] # user's two selected sequences
        st.write('You selected:')
        for record in sequences:
            if record.id in oligos_options:
                selected_sequences.append(record.seq)  # Append the selected sequences to the list 
                st.write(f">{record.id}")
                st.write(record.seq)
        
        sequence1 = Seq(selected_sequences[0])
        sequence2 = Seq(selected_sequences[1])
        
        # PCA:
        if option == 'Polymerase Cycling Assembly (PCA)':
            st.markdown('## <span style="font-size:30px;">Polymerase Cycling Assembly (PCA)</span>', unsafe_allow_html=True)
            output = pca_gibson.PCA(sequence1, sequence2, 20)
            if output == []:
                st.write('The selected sequences cannot be assembled using Polymerase Cycling Assembly (PCA)')
            else:
                st.write(output)
          
        # Restriction
        elif option == 'Restriction Assembly':
            st.markdown('## <span style="font-size:30px;">Restriction Assembly</span>', unsafe_allow_html=True)
            output = restriction_site.restriction(sequence1, sequence2)
            if output == []:
                st.write('The selected sequences cannot be assembled using Restriction assembly')
            else:
                st.write(output)
        
        # Gibson Assesmbly
        elif option == 'Gibson Assembly':
            st.markdown('## <span style="font-size:30px;">Gibson Assembly</span>', unsafe_allow_html=True)
            output = pca_gibson.Gibson(sequence1, sequence2, 15)
            if output == []:
                st.write('The selected sequences cannot be assembled using Gibson assembly')
            else:
                st.write(output)

        # Download the output as fasta file
        output_strings = []
        if output != []:
            if option == 'Polymerase Cycling Assembly (PCA)':
                for i, record in enumerate(output):
                    output_strings.append(f">pca{i + 1}\n{str(record)}")
            elif option == 'Restriction Assembly':
                for i, record in enumerate(output):
                    output_strings.append(f">restr{i + 1}\n{str(record)}")
            elif option == 'Gibson Assembly':
                for i, record in enumerate(output):
                    output_strings.append(f">gibson{i + 1}\n{str(record)}")
            
            output_content = "\n".join(output_strings)
            
            output_bytes = output_content.encode('utf-8')
            
            st.download_button(
                label="Download results as Fasta file",
                data=output_bytes,
                file_name='Results.fasta',
                mime='text/fasta',
            )
        
            
