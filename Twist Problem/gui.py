import PySimpleGUI as sg

# Define the background color (green pastel background)
background_color = '#D5E7B8' # green pastel background

font_size = 24

# Initial window layout of GUI
def create_initial_layout():
    return [
        [sg.Text("Please enter an oligo sequence:", background_color=background_color, text_color='black', font=('Any', font_size))],
        [sg.InputText(key='-SEQUENCE-', font=('Any', font_size))],  # Specify the font for the input bar
        [
            sg.Button("Next", button_color=('black', 'pink'), font=('Any', font_size)),
            sg.Button("Exit", button_color=('black', 'pink'), font=('Any', font_size))
        ]
    ]

# Convert oligo to DNA sequence
def oligo_to_dna(oligo):
    """
    Convert an oligonucleotide to a DNA sequence.

    Args:
        oligo (str): Oligonucleotide sequence using single-letter codes.

    Returns:
        str: DNA sequence.
    """
    conversion_dict = {
        'A': 'A',
        'C': 'C',
        'G': 'G',
        'T': 'T',
        'U': 'T',  # Assuming U (uracil) is converted to T (thymine)
        'R': 'AG',  # Purine (A or G)
        'Y': 'CT',  # Pyrimidine (C or T)
        'M': 'AC',  # Amino (A or C)
        'K': 'GT',  # Keto (G or T)
        'S': 'GC',  # Strong interaction (C or G)
        'W': 'AT',  # Weak interaction (A or T)
        'B': 'CGT',  # Not A (C or G or T)
        'D': 'AGT',  # Not C (A or G or T)
        'H': 'ACT',  # Not G (A or C or T)
        'V': 'ACG',  # Not T/U (A or C or G)
        'N': 'ACGT'  # Any base (A or C or G or T)
    }

    dna_sequence = ''.join(conversion_dict.get(base, '') for base in oligo.upper())
    return dna_sequence

# Create the initial window
window = sg.Window(title="Oligo App", layout=create_initial_layout(), margins=(300, 250), background_color=background_color)

show_result_layout = False  # Flag to determine whether to show the result layout

while True:
    event, values = window.read()
    
    # If user closes the window or clicks "Exit"
    if event == sg.WIN_CLOSED or event == "Exit":
        break
    
    # If user clicks Next
    if event == "Next" and not show_result_layout:
        oligo_sequence = values['-SEQUENCE-']
        
        # Convert the oligo sequence to DNA sequence (for now, just copying it)
        dna_sequence = oligo_to_dna(oligo_sequence)
        
        # Create a new layout to display the entered sequence
        result_layout = [
            [sg.Text(f"The converted DNA sequence is: {dna_sequence}", background_color=background_color, text_color='black', font=('Any', font_size))],
            [
                sg.Button("Try Another Oligo!", button_color=('black', 'pink'), font=('Any', font_size)),
                sg.Button("Exit", button_color=('black', 'pink'), font=('Any', font_size))
            ]
        ]
        
        # Close the current window
        window.close()
        
        # Create another window with the result layout
        window = sg.Window(title="Oligo App", layout=result_layout, margins=(300, 250), background_color=background_color)
        
        show_result_layout = True  # Set the flag to show result layout
    
    # If the user clicks 'Try another oligo'
    if event == "Try Another Oligo!" and show_result_layout:
        # Close the current window
        window.close()
        # Create a new instance of the initial layout
        window = sg.Window(title="Oligo App", layout=create_initial_layout(), margins=(300, 250), background_color=background_color)
        show_result_layout = False  # Reset the flag to show initial layout

window.close()
