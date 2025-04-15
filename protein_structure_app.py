import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio
import re
import pandas as pd
import numpy as np

# Page configuration
st.set_page_config(
    page_title="Protein Structure Predictor",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Apply custom CSS for improved styling
def apply_custom_css():
    st.markdown("""
    <style>
    .main .block-container {
        padding-top: 2rem;
        padding-bottom: 2rem;
    }
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
    }
    .stTabs [data-baseweb="tab"] {
        padding: 10px 16px;
        border-radius: 4px 4px 0px 0px;
    }
    .stTabs [aria-selected="true"] {
        background-color: #f0f0f0;
        font-weight: bold;
    }
    .stTextArea textarea {
        font-family: monospace;
    }
    .confidence-indicator {
        padding: 5px 10px;
        border-radius: 4px;
        font-weight: bold;
        display: inline-block;
    }
    h1 {
        margin-bottom: 0.5rem;
    }
    .stSidebar .block-container {
        padding-top: 2rem;
    }
    .info-box {
        background-color: #f8f9fa;
        border-left: 4px solid #4c78a8;
        padding: 1rem;
        margin-bottom: 1rem;
        border-radius: 4px;
    }
    .metric-card {
        background-color: #f8f9fa;
        padding: 1rem;
        border-radius: 8px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        text-align: center;
    }
    </style>
    """, unsafe_allow_html=True)

apply_custom_css()

# Function to render the 3D molecular structure
def render_mol(pdb, color_scheme='spectrum', spin=True):
    """
    Render the protein 3D structure using py3Dmol
    
    Parameters:
    -----------
    pdb : str
        PDB format data string
    color_scheme : str
        Color scheme for the 3D visualization
    spin : bool
        Whether to enable spinning animation
    """
    pdbview = py3Dmol.view(width=800, height=500)
    pdbview.addModel(pdb, 'pdb')
    
    # Apply different color schemes based on selection
    if color_scheme == "spectrum":
        pdbview.setStyle({'cartoon': {'color': 'spectrum'}})
    elif color_scheme == "chain":
        pdbview.setStyle({'cartoon': {'colorscheme': 'chainHetatm'}})
    elif color_scheme == "residue":
        pdbview.setStyle({'cartoon': {'colorscheme': 'amino'}})
    elif color_scheme == "secondary structure":
        pdbview.setStyle({'cartoon': {'color': 'secondary structure'}})
    elif color_scheme == "b-factor":
        pdbview.setStyle({'cartoon': {'colorbyfactor': 'b'}})
    
    pdbview.setBackgroundColor('white')
    pdbview.zoomTo()
    if spin:
        pdbview.spin(spin)
    showmol(pdbview, height=500, width=800)

# Function to evaluate the confidence level based on pLDDT score
def get_confidence(plddt):
    """
    Returns the confidence level based on pLDDT score
    
    Parameters:
    -----------
    plddt : float
        The pLDDT confidence score (0-100)
        
    Returns:
    --------
    str
        Confidence level description
    """
    if plddt >= 90:
        return "Very High"
    elif plddt >= 70:
        return "Confident"
    elif plddt >= 50:
        return "Low"
    else:
        return "Very Low"

# Function to get a color for the confidence level
def get_confidence_color(plddt):
    """
    Returns an appropriate color for the confidence level
    
    Parameters:
    -----------
    plddt : float
        The pLDDT confidence score (0-100)
        
    Returns:
    --------
    str
        Hex color code
    """
    if plddt >= 90:
        return "#198754"  # green
    elif plddt >= 70:
        return "#0d6efd"  # blue
    elif plddt >= 50:
        return "#fd7e14"  # orange
    else:
        return "#dc3545"  # red

# Function to fetch protein structure prediction from ESMFold API
def fetch_prediction(sequence):
    """
    Fetches protein structure prediction from ESMFold API
    
    Parameters:
    -----------
    sequence : str
        Amino acid sequence of the protein
        
    Returns:
    --------
    str
        PDB format data string
    """
    with st.spinner("🧬 Predicting protein structure... This may take a moment."):
        headers = {'Content-Type': 'application/x-www-form-urlencoded'}
        response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', 
                                headers=headers, data=sequence)
        if response.status_code == 200:
            return response.content.decode('utf-8')
        else:
            st.error(f"Error: Failed to predict structure. Status code: {response.status_code}")
            return None

# App header with logo and title
st.markdown("""
<div style="display: flex; align-items: center; margin-bottom: 1rem;">
    <div style="font-size: 3rem; margin-right: 1rem;">🧬</div>
    <div>
        <h1 style="margin: 0; padding: 0;">Protein Structure Predictor</h1>
        <p style="margin: 0; color: #666;">Powered by ESMFold</p>
    </div>
</div>
""", unsafe_allow_html=True)

st.markdown('<hr style="margin: 1rem 0">', unsafe_allow_html=True)

# Sidebar with input options
with st.sidebar:
    st.markdown("""
    <div style="text-align: center; margin-bottom: 1rem;">
        <div style="font-size: 2rem;">🔬</div>
        <h3 style="margin-top: 0.5rem;">Protein Input & Settings</h3>
    </div>
    """, unsafe_allow_html=True)
    
    # Example sequences
    example_sequences = {
        "Current sequence": "",
        "T4 Lysozyme (small)": "MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL",
        "GFP (medium)": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
        "Human Hemoglobin (complex)": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
    }

    selected_example = st.selectbox(
        "Select an example protein",
        options=list(example_sequences.keys()),
        help="Choose from pre-defined protein sequences"
    )
    
    sequence_input = example_sequences[selected_example]
    
    if selected_example == "Current sequence":
        sequence_input = st.session_state.get("sequence", "")
        st.session_state.sequence = sequence_input
    else:
        st.session_state.sequence = sequence_input
    
    sequence = st.text_area(
        "Enter Protein Sequence",
        sequence_input,
        height=250,
        help="Paste your amino acid sequence here. Use standard one-letter amino acid codes."
    )
    
    st.markdown('<hr style="margin: 1rem 0">', unsafe_allow_html=True)
    
    # Visualization options
    st.markdown("""
    <div style="text-align: center; margin-bottom: 1rem;">
        <div style="font-size: 2rem;">🔧</div>
        <h3 style="margin-top: 0.5rem;">Visualization Options</h3>
    </div>
    """, unsafe_allow_html=True)
    
    color_scheme = st.selectbox(
        "Color Scheme",
        ["spectrum", "chain", "residue", "secondary structure", "b-factor"],
        help="Choose different coloring methods for the 3D structure"
    )
    
    spin = st.checkbox(
        "Spin Structure",
        value=True,
        help="Toggle rotation animation for the 3D model"
    )
    
    # Action buttons
    predict = st.button("🔍 Predict Structure", type="primary", use_container_width=True)
    reset = st.button("🔄 Reset", type="secondary", use_container_width=True)
    
    if reset:
        for key in list(st.session_state.keys()):
            del st.session_state[key]
        st.rerun()
    
    # About box
    st.markdown("""
    <div class="info-box">
        <h4 style="margin-top: 0;">About</h4>
        <p>This app uses ESMFold to predict protein structures from amino acid sequences. 
        The model produces a confidence score (pLDDT) for each prediction.</p>
    </div>
    """, unsafe_allow_html=True)

# Main content area
if predict and sequence:
    # Clean the sequence to remove whitespace and invalid characters
    clean_sequence = re.sub(r'[^A-Za-z]', '', sequence)
    
    if not clean_sequence:
        st.error("Please enter a valid amino acid sequence.")
    else:
        pdb_string = fetch_prediction(clean_sequence)
        
        if pdb_string:
            # Save PDB file
            with open('predicted.pdb', 'w') as f:
                f.write(pdb_string)
            
            # Load structure and calculate metrics
            struct = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
            plddt = round(struct.b_factor.mean(), 2)
            seq_length = len(clean_sequence)
            confidence = get_confidence(plddt)
            confidence_color = get_confidence_color(plddt)
            
            # Store pLDDT values in session state
            st.session_state.plddts = struct.b_factor.tolist()
            
            # Create tabs for different views
            tab1, tab2, tab3 = st.tabs([
                "📊 3D Structure", 
                "🧪 Sequence Analysis", 
                "📄 PDB Data"
            ])
            
            # Tab 1: 3D Structure
            with tab1:
                st.markdown("""
                <h3 style="margin-bottom: 1.5rem;">Predicted Protein Structure</h3>
                """, unsafe_allow_html=True)
                
                # Show the 3D structure
                render_mol(pdb_string, color_scheme=color_scheme, spin=spin)
                
                # Metrics in cards
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.markdown(f"""
                    <div class="metric-card">
                        <div style="font-size: 0.8rem; color: #666;">Sequence Length</div>
                        <div style="font-size: 1.8rem; font-weight: bold;">{seq_length}</div>
                        <div style="font-size: 0.8rem; color: #666;">amino acids</div>
                    </div>
                    """, unsafe_allow_html=True)
                
                with col2:
                    st.markdown(f"""
                    <div class="metric-card">
                        <div style="font-size: 0.8rem; color: #666;">Average pLDDT</div>
                        <div style="font-size: 1.8rem; font-weight: bold;">{plddt}</div>
                        <div style="font-size: 0.8rem; color: #666;">confidence score</div>
                    </div>
                    """, unsafe_allow_html=True)
                
                with col3:
                    st.markdown(f"""
                    <div class="metric-card">
                        <div style="font-size: 0.8rem; color: #666;">Confidence</div>
                        <div style="font-size: 1.8rem; font-weight: bold; color: {confidence_color};">{confidence}</div>
                        <div style="font-size: 0.8rem; color: #666;">prediction reliability</div>
                    </div>
                    """, unsafe_allow_html=True)
                
                # Visualization explanation
                st.markdown("""
                <div class="info-box" style="margin-top: 1.5rem;">
                    <h4 style="margin-top: 0;">Structure Visualization</h4>
                    <p>
                        The 3D model shows the predicted protein structure. You can:
                        <ul>
                            <li>Rotate: Click and drag</li>
                            <li>Zoom: Scroll wheel</li>
                            <li>Change colors: Use the sidebar options</li>
                        </ul>
                    </p>
                </div>
                """, unsafe_allow_html=True)
            
            # Tab 2: Sequence Analysis
            with tab2:
                st.markdown("""
                <h3 style="margin-bottom: 1.5rem;">Sequence Analysis</h3>
                """, unsafe_allow_html=True)
                
                # pLDDT confidence visualization
                if st.session_state.plddts:
                    st.markdown("""
                    <h4>pLDDT Confidence Scores</h4>
                    <p>
                        pLDDT (predicted Local Distance Difference Test) is a per-residue confidence score.
                        Higher values indicate more reliable predictions:
                    </p>
                    """, unsafe_allow_html=True)
                    
                    # Confidence legend
                    col1, col2, col3, col4 = st.columns(4)
                    col1.markdown(f'<div class="confidence-indicator" style="background-color: {get_confidence_color(95)}; color: white;">Very High (90+)</div>', unsafe_allow_html=True)
                    col2.markdown(f'<div class="confidence-indicator" style="background-color: {get_confidence_color(80)}; color: white;">High (70-90)</div>', unsafe_allow_html=True)
                    col3.markdown(f'<div class="confidence-indicator" style="background-color: {get_confidence_color(60)}; color: white;">Medium (50-70)</div>', unsafe_allow_html=True)
                    col4.markdown(f'<div class="confidence-indicator" style="background-color: {get_confidence_color(40)}; color: white;">Low (<50)</div>', unsafe_allow_html=True)
                    
                    # Create pLDDT chart
                    plddt_values = st.session_state.plddts
                    residue_positions = list(range(1, len(plddt_values) + 1))
                    
                    # Create enhanced chart with color bands
                    plddt_chart_data = pd.DataFrame({
                        "residue": residue_positions,
                        "pLDDT": plddt_values
                    })
                    
                    # Create a more visually appealing chart
                    st.line_chart(
                        plddt_chart_data.set_index("residue"),
                        use_container_width=True,
                        height=300
                    )
                
                # Amino Acid Composition
                st.markdown("""
                <h4>Amino Acid Composition</h4>
                <p>Distribution of amino acids in the sequence:</p>
                """, unsafe_allow_html=True)
                
                if clean_sequence:
                    aa_counts = {}
                    for aa in clean_sequence:
                        aa_counts[aa] = aa_counts.get(aa, 0) + 1
                    
                    # Sort amino acids by chemical properties for better visualization
                    aa_groups = {
                        'Hydrophobic': ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'],
                        'Polar': ['S', 'T', 'N', 'Q'],
                        'Positively Charged': ['K', 'R', 'H'],
                        'Negatively Charged': ['D', 'E'],
                        'Special': ['C', 'G', 'P']
                    }
                    
                    # Convert to dataframe
                    aa_data = pd.DataFrame({
                        "Amino Acid": list(aa_counts.keys()),
                        "Count": list(aa_counts.values())
                    })
                    
                    # Add percentage
                    aa_data['Percentage'] = (aa_data['Count'] / aa_data['Count'].sum() * 100).round(1)
                    
                    # Sort by count
                    aa_data = aa_data.sort_values('Count', ascending=False)
                    
                    # Display bar chart
                    st.bar_chart(aa_data.set_index("Amino Acid")['Count'], use_container_width=True)
                    
                    # Display table with counts and percentages
                    st.markdown("""
                    <h4>Amino Acid Details</h4>
                    """, unsafe_allow_html=True)
                    
                    aa_data['Percentage'] = aa_data['Percentage'].astype(str) + '%'
                    st.dataframe(aa_data, use_container_width=True)
                    
                    # Hydrophobicity analysis
                    hydrophobic_count = sum([aa_counts.get(aa, 0) for aa in aa_groups['Hydrophobic']])
                    polar_count = sum([aa_counts.get(aa, 0) for aa in aa_groups['Polar']])
                    positive_count = sum([aa_counts.get(aa, 0) for aa in aa_groups['Positively Charged']])
                    negative_count = sum([aa_counts.get(aa, 0) for aa in aa_groups['Negatively Charged']])
                    special_count = sum([aa_counts.get(aa, 0) for aa in aa_groups['Special']])
                    
                    total = len(clean_sequence)
                    
                    # Create properties chart
                    prop_data = pd.DataFrame({
                        'Property': ['Hydrophobic', 'Polar', 'Positively Charged', 'Negatively Charged', 'Special'],
                        'Count': [hydrophobic_count, polar_count, positive_count, negative_count, special_count],
                        'Percentage': [
                            round(hydrophobic_count/total*100, 1),
                            round(polar_count/total*100, 1),
                            round(positive_count/total*100, 1),
                            round(negative_count/total*100, 1),
                            round(special_count/total*100, 1)
                        ]
                    })
                    
                    # Show property distribution
                    st.markdown("""
                    <h4>Amino Acid Property Distribution</h4>
                    <p>Amino acids grouped by their chemical properties:</p>
                    """, unsafe_allow_html=True)
                    
                    fig_data = prop_data.set_index('Property')
                    st.bar_chart(fig_data['Percentage'], use_container_width=True)
                    
                    # Property descriptions
                    st.markdown("""
                    <div class="info-box">
                        <h4 style="margin-top: 0;">Property Descriptions</h4>
                        <ul>
                            <li><strong>Hydrophobic:</strong> Non-polar amino acids that tend to cluster in the protein core (A, V, I, L, M, F, Y, W)</li>
                            <li><strong>Polar:</strong> Polar, uncharged amino acids that interact with water (S, T, N, Q)</li>
                            <li><strong>Positively Charged:</strong> Basic amino acids with a positive charge (K, R, H)</li>
                            <li><strong>Negatively Charged:</strong> Acidic amino acids with a negative charge (D, E)</li>
                            <li><strong>Special:</strong> Amino acids with special structural roles (C, G, P)</li>
                        </ul>
                    </div>
                    """, unsafe_allow_html=True)
            
            # Tab 3: PDB Data
            with tab3:
                st.markdown("""
                <h3 style="margin-bottom: 1.5rem;">PDB Format Data</h3>
                <p>The Protein Data Bank (PDB) format contains atomic coordinates and other information about the protein structure.</p>
                """, unsafe_allow_html=True)
                
                # Download button with improved styling
                st.markdown("""
                <div style="margin-bottom: 1rem;">
                    <p>Download the PDB file to use with other molecular visualization software:</p>
                </div>
                """, unsafe_allow_html=True)
                
                st.download_button(
                    "💾 Download PDB File",
                    pdb_string,
                    file_name="predicted.pdb",
                    mime="text/plain",
                    use_container_width=True
                )
                
                # Display PDB data in a formatted text area
                st.markdown("""
                <h4 style="margin-top: 1.5rem;">PDB File Content</h4>
                <p>This is the raw PDB format data:</p>
                """, unsafe_allow_html=True)
                
                st.text_area(
                    "",
                    pdb_string,
                    height=300,
                    help="Copy of the PDB file content"
                )
                
                # PDB format explanation
                st.markdown("""
                <div class="info-box">
                    <h4 style="margin-top: 0;">About PDB Format</h4>
                    <p>
                        The PDB format contains:
                        <ul>
                            <li>Atomic coordinates for all atoms in the structure</li>
                            <li>Secondary structure annotations</li>
                            <li>B-factor values (pLDDT scores in this case)</li>
                            <li>Connection information between atoms</li>
                        </ul>
                    </p>
                </div>
                """, unsafe_allow_html=True)

else:
    # Welcome state
    st.markdown("""
    <div style="text-align: center; padding: 3rem 1rem;">
        <img src="https://cdn-icons-png.flaticon.com/512/2942/2942789.png" style="width: 120px; margin-bottom: 1.5rem;" />
        <h2>Welcome to the Protein Structure Predictor</h2>
        <p style="max-width: 600px; margin: 1rem auto; font-size: 1.1rem; color: #555;">
            This tool uses AI to predict the 3D structure of proteins from their amino acid sequences.
            Select an example protein or enter your own sequence to get started.
        </p>
        
        <div style="background-color: #f8f9fa; border-radius: 8px; padding: 1.5rem; max-width: 600px; margin: 0 auto 1.5rem auto; text-align: left;">
            <h3 style="margin-top: 0;">About Protein Structure</h3>
            <p>
                Proteins are complex molecules that play critical roles in all biological processes. Their function is 
                determined by their three-dimensional structure, which is formed by folding chains of amino acids.
            </p>
            <p>
                <strong>Key structural levels:</strong>
            </p>
            <ul>
                <li><strong>Primary structure:</strong> The sequence of amino acids</li>
                <li><strong>Secondary structure:</strong> Local folded structures (α-helices and β-sheets)</li>
                <li><strong>Tertiary structure:</strong> The overall 3D arrangement of the entire protein</li>
                <li><strong>Quaternary structure:</strong> The arrangement of multiple protein subunits</li>
            </ul>
            <p>
                Understanding protein structure helps scientists develop new drugs, engineer proteins with enhanced functions,
                and understand disease mechanisms at the molecular level.
            </p>
        </div>
        
        <div style="background-color: #f8f9fa; border-radius: 8px; padding: 1.5rem; max-width: 600px; margin: 0 auto; text-align: left;">
            <h3 style="margin-top: 0;">How to use:</h3>
            <ol style="margin-bottom: 0;">
                <li>Select an example protein or paste your own sequence in the sidebar</li>
                <li>Click the "Predict Structure" button</li>
                <li>Explore the 3D structure, sequence analysis, and download the PDB file</li>
            </ol>
        </div>
    </div>
    """, unsafe_allow_html=True)
