import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio
import re
import pandas as pd
import plotly.express as px
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import time
import os
from io import StringIO

# ====== CUSTOM CSS & PAGE CONFIG ======
st.set_page_config(
    layout='wide',
    page_title="ESMFold Protein Predictor",
    page_icon="🧬",
    initial_sidebar_state="expanded"
)

# Custom CSS for modern UI
st.markdown("""
<style>
    /* Main app background */
    .stApp {
        background-color: #f8fafc;
    }
    
    /* Sidebar styling */
    [data-testid="stSidebar"] {
        background: linear-gradient(180deg, #2563eb 0%, #1e40af 100%) !important;
        padding: 1.5rem !important;
    }
    
    /* Sidebar text colors */
    .sidebar .sidebar-content {
        color: white !important;
    }
    
    .sidebar .stSelectbox label, 
    .sidebar .stCheckbox label,
    .sidebar .stTextArea label,
    .sidebar .stMarkdown h1 {
        color: white !important;
    }
    
    /* Button styling */
    div.stButton > button:first-child {
        background-color: #2563eb;
        color: white;
        border-radius: 8px;
        padding: 0.5rem 1rem;
        border: none;
        transition: all 0.3s;
        font-weight: 500;
    }
    
    div.stButton > button:first-child:hover {
        background-color: #1e40af !important;
        transform: translateY(-1px);
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
    }
    
    /* Text area styling */
    .stTextArea textarea {
        border-radius: 12px !important;
        padding: 1rem !important;
        border: 1px solid #e2e8f0 !important;
    }
    
    /* Tabs styling */
    .stTabs [role="tablist"] {
        gap: 8px;
    }
    
    .stTabs [role="tab"] {
        border-radius: 8px 8px 0 0 !important;
        padding: 0.5rem 1rem !important;
        background: #e2e8f0 !important;
        transition: all 0.3s;
    }
    
    .stTabs [role="tab"][aria-selected="true"] {
        background: #2563eb !important;
        color: white !important;
    }
    
    /* Metric cards */
    [data-testid="stMetric"] {
        background: white;
        border-radius: 12px;
        padding: 1rem;
        box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        border: 1px solid #e2e8f0;
    }
    
    /* Custom header */
    .custom-header {
        background: linear-gradient(90deg, #2563eb 0%, #1e40af 100%);
        padding: 2rem;
        border-radius: 12px;
        color: white;
        margin-bottom: 1.5rem;
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
    }
    
    /* Footer styling */
    .custom-footer {
        text-align: center;
        padding: 1.5rem;
        margin-top: 3rem;
        color: #64748b;
        font-size: 0.9rem;
    }
    
    /* Progress bar styling */
    .stProgress > div > div > div > div {
        background-color: #2563eb;
    }
    
    /* Tooltip styling */
    .stTooltip {
        font-size: 0.9rem !important;
    }
</style>
""", unsafe_allow_html=True)

# ====== HEADER SECTION ======
st.markdown("""
<div class="custom-header">
    <h1 style="margin: 0; color: white;">🧬 Protein Structure Predictor</h1>
    <p style="margin: 0.5rem 0 0; font-size: 1.1rem; opacity: 0.9;">
        Predict 3D protein structures using ESMFold's cutting-edge AI
    </p>
</div>
""", unsafe_allow_html=True)

# ====== ENHANCED FUNCTIONS ======
def render_mol(pdb, color_scheme='spectrum', spin=True, style='cartoon', width=800, height=500):
    pdbview = py3Dmol.view(width=width, height=height)
    pdbview.addModel(pdb, 'pdb')
    
    if style == 'cartoon':
        pdbview.setStyle({'cartoon': {'color': color_scheme}})
    elif style == 'stick':
        pdbview.setStyle({'stick': {'color': color_scheme}})
    elif style == 'sphere':
        pdbview.setStyle({'sphere': {'color': color_scheme}})
    elif style == 'surface':
        pdbview.addSurface(py3Dmol.VDW, {'opacity': 0.7, 'color': color_scheme})
        pdbview.setStyle({'cartoon': {'color': 'spectrum'}})
    
    pdbview.setBackgroundColor('white')
    pdbview.zoomTo()
    pdbview.spin(spin)
    showmol(pdbview, height=height, width=width)

def get_confidence(plddt):
    if plddt >= 90:
        return "Very High", "#10B981"
    elif plddt >= 70:
        return "Confident", "#3B82F6"
    elif plddt >= 50:
        return "Low", "#F59E0B"
    else:
        return "Very Low", "#EF4444"

def fetch_prediction(sequence):
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence)
    return response.content.decode('utf-8')

def analyze_protein(sequence):
    """Perform advanced protein sequence analysis"""
    analysis = ProteinAnalysis(sequence)
    
    # Calculate basic properties
    molecular_weight = analysis.molecular_weight()
    aromaticity = analysis.aromaticity()
    instability_index = analysis.instability_index()
    isoelectric_point = analysis.isoelectric_point()
    secondary_structure = analysis.secondary_structure_fraction()
    
    # Calculate amino acid percentages
    aa_percent = analysis.get_amino_acids_percent()
    
    return {
        'molecular_weight': molecular_weight,
        'aromaticity': aromaticity,
        'instability_index': instability_index,
        'isoelectric_point': isoelectric_point,
        'secondary_structure': secondary_structure,
        'aa_percent': aa_percent
    }

def validate_sequence(sequence):
    """Validate protein sequence input"""
    sequence = re.sub(r'\s+', '', sequence.upper())
    if not sequence:
        return False, "Sequence is empty"
    
    valid_aas = set('ACDEFGHIKLMNPQRSTVWY')
    invalid_chars = set(sequence) - valid_aas
    
    if invalid_chars:
        return False, f"Invalid amino acids: {', '.join(invalid_chars)}"
    
    return True, sequence

# ====== SIDEBAR CONTROLS ======
with st.sidebar:
    # Sidebar header with logo placeholder
    st.markdown("""
    <div style="text-align: center; margin-bottom: 2rem;">
        <div style="font-size: 2rem; margin-bottom: 0.5rem;">🧪</div>
        <h2 style="color: white; margin: 0;">Protein Input</h2>
    </div>
    """, unsafe_allow_html=True)
    
    example_sequences = {
        "Select example": "",
        "T4 Lysozyme (small)": "MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL",
        "GFP (medium)": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
        "Human Hemoglobin (complex)": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
        "Insulin": "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
    }

    selected_example = st.selectbox(
        "Select an example protein", 
        options=list(example_sequences.keys()),
        key="protein_select"
    )
    
    sequence_input = example_sequences[selected_example]
    
    if selected_example == "Select example":
        sequence_input = st.session_state.get("sequence", "")
        st.session_state.sequence = sequence_input
    else:
        st.session_state.sequence = sequence_input

    sequence = st.text_area(
        "Enter Protein Sequence", 
        sequence_input, 
        height=250,
        key="seq_input",
        help="Enter a protein sequence using standard amino acid codes (A-Z)"
    )

    st.markdown("---")
    st.title("🔧 Visualization Options")
    
    col1, col2 = st.columns(2)
    with col1:
        color_scheme = st.selectbox(
            "Color Scheme", 
            ["spectrum", "chain", "residue", "secondary structure"],
            key="color_scheme"
        )
        
        style = st.selectbox(
            "Visualization Style",
            ["cartoon", "stick", "sphere", "surface"],
            key="style"
        )
    
    with col2:
        spin = st.checkbox("Spin Structure", value=True, key="spin")
        show_sidechains = st.checkbox("Show Side Chains", value=False, key="sidechains")
    
    st.markdown("---")
    st.title("⚙️ Advanced Options")
    max_retries = st.slider("Max API Retries", 1, 5, 3, help="Number of retries if API fails")
    timeout = st.slider("Timeout (seconds)", 10, 120, 30, help="API request timeout")
    
    reset = st.button("Reset", type="primary", key="reset")
    predict = st.button("Predict Structure", key="predict")
    
    if reset:
        for key in st.session_state.keys():
            del st.session_state[key]
        st.rerun()

# ====== MAIN CONTENT AREA ======
if predict:
    # Validate sequence
    is_valid, validation_msg = validate_sequence(sequence)
    
    if not is_valid:
        st.error(f"Invalid sequence: {validation_msg}")
        st.stop()
    
    clean_sequence = validation_msg  # This is the validated sequence
    
    if len(clean_sequence) > 400:
        st.warning("⚠️ Sequences longer than 400 amino acids may take longer to process and have lower accuracy.")
    
    with st.spinner('🚀 Predicting protein structure...'):
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        pdb_string = None
        retry_count = 0
        
        while retry_count < max_retries and not pdb_string:
            try:
                status_text.text(f"Attempt {retry_count + 1} of {max_retries}...")
                progress_bar.progress((retry_count + 1) / max_retries)
                
                start_time = time.time()
                pdb_string = fetch_prediction(clean_sequence)
                elapsed_time = time.time() - start_time
                
                # Save to session state
                st.session_state.pdb_string = pdb_string
                st.session_state.sequence = clean_sequence
                
                # Save to file
                with open('predicted.pdb', 'w') as f:
                    f.write(pdb_string)
                
                # Parse structure
                struct = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
                plddt = round(struct.b_factor.mean(), 2)
                seq_length = len(clean_sequence)
                confidence, confidence_color = get_confidence(plddt)
                
                st.session_state.plddts = struct.b_factor.tolist()
                st.session_state.struct = struct
                st.session_state.analysis = analyze_protein(clean_sequence)
                
                break
                
            except Exception as e:
                retry_count += 1
                if retry_count == max_retries:
                    st.error(f"Failed to predict structure after {max_retries} attempts. Error: {str(e)}")
                    st.stop()
                time.sleep(2)  # Wait before retrying

        progress_bar.empty()
        status_text.empty()
        
        # Show success message
        st.success(f"Structure predicted successfully in {elapsed_time:.2f} seconds!")
        
        tab1, tab2, tab3, tab4 = st.tabs(["3D Structure", "Sequence Analysis", "PDB Data", "Advanced Analysis"])

        with tab1:
            col1, col2 = st.columns([3, 1])
            
            with col1:
                st.subheader("Predicted Protein Structure")
                render_mol(
                    pdb_string, 
                    color_scheme=color_scheme, 
                    spin=spin,
                    style=style,
                    width=800,
                    height=600
                )
            
            with col2:
                st.subheader("Prediction Metrics")
                
                # Create a metric for each important value
                st.metric("Sequence Length", seq_length)
                
                # pLDDT with colored indicator
                confidence_text = f"<span style='color:{confidence_color}; font-weight:bold;'>{confidence}</span>"
                st.metric(
                    "Average pLDDT", 
                    plddt,
                    help="pLDDT score indicates prediction confidence (0-100)"
                )
                st.markdown(f"**Confidence:** {confidence_text}", unsafe_allow_html=True)
                
                # Molecular weight if available
                if 'analysis' in st.session_state:
                    st.metric(
                        "Molecular Weight", 
                        f"{st.session_state.analysis['molecular_weight']/1000:.2f} kDa"
                    )
                
                # Download buttons
                st.download_button(
                    "Download PDB", 
                    pdb_string, 
                    file_name="predicted_structure.pdb", 
                    mime="text/plain"
                )
                
                st.download_button(
                    "Download pLDDT Data",
                    pd.DataFrame({'Residue': range(1, len(st.session_state.plddts)+1),
                                'pLDDT': st.session_state.plddts}).to_csv(index=False),
                    file_name="plddt_scores.csv",
                    mime="text/csv"
                )

        with tab2:
            st.subheader('Sequence Analysis')
            
            if st.session_state.plddts:
                # Enhanced pLDDT plot with Plotly
                st.markdown("""
                ### pLDDT per Residue
                pLDDT (predicted Local Distance Difference Test) estimates confidence per residue:
                - 🟢 **90+** : Very high confidence
                - 🟡 **70-90**: High confidence
                - 🟠 **50-70**: Medium confidence
                - 🔴 **Below 50**: Low confidence
                """)
                
                plddt_data = pd.DataFrame({
                    "Residue": list(range(1, len(st.session_state.plddts) + 1)), 
                    "pLDDT": st.session_state.plddts,
                    "Confidence": [get_confidence(score)[0] for score in st.session_state.plddts]
                })
                
                fig = px.line(
                    plddt_data, 
                    x="Residue", 
                    y="pLDDT",
                    color="Confidence",
                    color_discrete_map={
                        "Very High": "#10B981",
                        "Confident": "#3B82F6",
                        "Low": "#F59E0B",
                        "Very Low": "#EF4444"
                    },
                    title="Residue-wise Confidence Scores",
                    labels={"pLDDT": "Confidence Score"},
                    height=400
                )
                
                # Add horizontal lines for confidence thresholds
                fig.add_hline(y=90, line_dash="dot", line_color="#10B981", opacity=0.3)
                fig.add_hline(y=70, line_dash="dot", line_color="#3B82F6", opacity=0.3)
                fig.add_hline(y=50, line_dash="dot", line_color="#F59E0B", opacity=0.3)
                
                st.plotly_chart(fig, use_container_width=True)
            
            # Amino acid composition with enhanced visualization
            if clean_sequence:
                clean_seq = re.sub(r'[^A-Za-z]', '', clean_sequence)
                aa_counts = {}
                for aa in clean_seq:
                    aa_counts[aa] = aa_counts.get(aa, 0) + 1
                
                st.markdown("### Amino Acid Composition")
                
                # Create two columns for the charts
                col1, col2 = st.columns(2)
                
                with col1:
                    # Bar chart
                    aa_data = pd.DataFrame({
                        "Amino Acid": list(aa_counts.keys()), 
                        "Count": list(aa_counts.values())
                    }).sort_values("Count", ascending=False)
                    
                    fig = px.bar(
                        aa_data, 
                        x="Amino Acid", 
                        y="Count",
                        title="Amino Acid Counts",
                        color="Amino Acid",
                        color_discrete_sequence=px.colors.qualitative.Plotly
                    )
                    st.plotly_chart(fig, use_container_width=True)
                
                with col2:
                    # Pie chart
                    fig = px.pie(
                        aa_data,
                        names="Amino Acid",
                        values="Count",
                        title="Amino Acid Distribution",
                        hole=0.3
                    )
                    st.plotly_chart(fig, use_container_width=True)

        with tab3:
            st.subheader("PDB Format Data")
            
            col1, col2 = st.columns([3, 1])
            with col1:
                st.download_button(
                    "Download PDB File", 
                    pdb_string, 
                    file_name="predicted.pdb", 
                    mime="text/plain"
                )
            
            with col2:
                st.download_button(
                    "Download FASTA", 
                    f">predicted_protein\n{clean_sequence}", 
                    file_name="sequence.fasta", 
                    mime="text/plain"
                )
            
            st.code(pdb_string, language='pdb')

        with tab4:
            st.subheader("Advanced Protein Analysis")
            
            if 'analysis' in st.session_state:
                analysis = st.session_state.analysis
                
                st.markdown("### Protein Properties")
                cols = st.columns(4)
                with cols[0]:
                    st.metric("Molecular Weight", f"{analysis['molecular_weight']/1000:.2f} kDa")
                with cols[1]:
                    st.metric("Aromaticity", f"{analysis['aromaticity']:.3f}")
                with cols[2]:
                    st.metric("Instability Index", f"{analysis['instability_index']:.2f}")
                with cols[3]:
                    st.metric("Isoelectric Point", f"{analysis['isoelectric_point']:.2f}")
                
                st.markdown("### Secondary Structure Prediction")
                ss_data = pd.DataFrame({
                    "Type": ["Helix", "Turn", "Sheet"],
                    "Fraction": [
                        analysis['secondary_structure'][0],
                        analysis['secondary_structure'][1],
                        analysis['secondary_structure'][2]
                })
                
                fig = px.bar(
                    ss_data,
                    x="Type",
                    y="Fraction",
                    title="Secondary Structure Fractions",
                    color="Type",
                    color_discrete_sequence=px.colors.qualitative.Pastel
                )
                st.plotly_chart(fig, use_container_width=True)
                
                st.markdown("### Amino Acid Percentages")
                aa_percent = pd.DataFrame({
                    "Amino Acid": list(analysis['aa_percent'].keys()),
                    "Percentage": list(analysis['aa_percent'].values())
                }).sort_values("Percentage", ascending=False)
                
                fig = px.bar(
                    aa_percent,
                    x="Amino Acid",
                    y="Percentage",
                    title="Amino Acid Percentages",
                    color="Amino Acid",
                    color_discrete_sequence=px.colors.qualitative.Plotly
                )
                st.plotly_chart(fig, use_container_width=True)

elif not predict and not st.session_state.get("sequence"):
    st.info("👈 Enter a protein sequence and click **Predict Structure** to begin")

# ====== FOOTER ======
st.markdown("""
<div class="custom-footer">
    <hr style="border: 0.5px solid #e2e8f0; margin: 1.5rem 0;">
    <p>Protein Structure Predictor v2.0 | Powered by ESMFold API</p>
    <p style="font-size: 0.8rem; margin-top: 0.5rem;">
        For research use only | Not for clinical or diagnostic use
    </p>
</div>
""", unsafe_allow_html=True)
