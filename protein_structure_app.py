import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio
import re
import pandas as pd

# ====== SESSION STATE INITIALIZATION ======
for key, default_value in {
    "prediction_made": False,
    "sequence": "",
    "color_scheme": "spectrum",
    "spin": True,
    "plddts": []
}.items():
    if key not in st.session_state:
        st.session_state[key] = default_value

# ====== CUSTOM CSS & PAGE CONFIG ======
st.set_page_config(
    layout='wide',
    page_title="ESMFold Protein Predictor",
    page_icon="🧬",
    initial_sidebar_state="expanded"
)

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;700&display=swap');

body, .stApp, [data-testid="stSidebar"] {
    font-family: 'Inter', sans-serif !important;
}

/* Main app background */
.stApp {
    background-color: #f8fafc;
}

/* Sidebar styling */
[data-testid="stSidebar"] {
    background: linear-gradient(180deg, #2563eb 0%, #1e40af 100%) !important;
    padding: 1.1rem 1.1rem 1.1rem 1.1rem !important;
    color: #f0f0f0 !important;
    min-width: 340px !important;
    max-width: 340px !important;
    position: relative;
}

/* Remove the white top bar at the top of the sidebar for seamless join */
[data-testid="stSidebar"]::before {
    background-color: transparent !important;
    height: 0 !important;
}

/* Adjust sidebar content to not overlap the top */
[data-testid="stSidebar"] > div:first-child {
    margin-top: 0 !important;
}

/* Sidebar widget labels and inputs */
[data-testid="stSidebar"] label, 
[data-testid="stSidebar"] .stTextArea textarea, 
[data-testid="stSidebar"] .stTextInput input, 
[data-testid="stSidebar"] select,
[data-testid="stSidebar"] .stCheckbox label {
    color: #fff !important;
    font-weight: 600;
}

/* Make the Spin Structure label bold and white */
[data-testid="stSidebar"] .stCheckbox label {
    color: #fff !important;
    font-weight: 700 !important;
}

/* Sidebar input fields background and border */
[data-testid="stSidebar"] .stTextArea textarea, 
[data-testid="stSidebar"] .stTextInput input, 
[data-testid="stSidebar"] select {
    background-color: #1e3a8a !important;
    border: 1px solid #3b82f6 !important;
    color: white !important;
    border-radius: 8px;
    padding: 0.5rem;
    font-size: 1rem;
}

/* Sidebar buttons */
[data-testid="stSidebar"] div.stButton > button {
    background-color: #2563eb !important;
    color: white !important;
    border-radius: 8px !important;
    padding: 0.5rem 1rem !important;
    font-weight: 600 !important;
    transition: background-color 0.3s ease;
    border: none !important;
    cursor: pointer;
    width: 100%;
    margin-top: 0.25rem;
    margin-bottom: 0.25rem;
}
[data-testid="stSidebar"] div.stButton > button:hover {
    background-color: #1e40af !important;
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
    font-weight: 600;
    font-family: 'Inter', sans-serif !important;
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
    font-family: 'Inter', sans-serif !important;
}

/* Text area styling */
.stTextArea textarea {
    border-radius: 12px !important;
    padding: 0.65rem !important;
    border: 1px solid #e2e8f0 !important;
    font-family: 'Inter', sans-serif !important;
    min-height: 110px !important;
    max-height: 140px !important;
    font-size: 1.05rem !important;
    margin-bottom: 0.3rem !important;
}

/* Custom header - joined with sidebar, no left/top radius, lifted up */
.custom-header {
    background: linear-gradient(90deg, #2563eb 0%, #1e40af 100%);
    margin-top: 0 !important;
    margin-left: 0 !important;
    padding-top: 0.7rem !important;
    padding-bottom: 0.7rem !important;
    padding-left: 1.5rem !important;
    padding-right: 1.5rem !important;
    border-radius: 0 0 12px 12px !important; /* No top or left rounding */
    color: white;
    margin-bottom: 0.5rem;
    box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.08);
    font-family: 'Inter', sans-serif !important;
    z-index: 10;
    position: relative;
}

/* Remove default Streamlit main block top padding */
.block-container {
    padding-top: 0rem !important;
}

/* Reduce space below the protein cards and above the footer */
.main-content {
    position: relative;
    z-index: 1;
    margin-bottom: 0.5rem !important;
}
.protein-card {
    margin-bottom: 0.5rem !important;
}
.custom-footer {
    margin-top: 0.5rem !important;
    padding: 1rem 0 !important;
    color: #64748b;
    font-size: 0.9rem;
    font-family: 'Inter', sans-serif !important;
    text-align: center;
}
hr {margin: 0.7rem 0;}
/* Protein card and other unchanged styles omitted for brevity, keep your existing ones */
</style>
""", unsafe_allow_html=True)

# ====== PROTEIN TEMPLATE DATA ======
protein_templates = {
    "T4 Lysozyme (small)": {
        "sequence": "MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL",
        "image_url": "https://raw.githubusercontent.com/2-0-0-5-emil/protein-structure-app/main/lyzo.jpeg"
    },
    "GFP (medium)": {
        "sequence": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
        "image_url": "https://raw.githubusercontent.com/2-0-0-5-emil/protein-structure-app/main/gfp.jpeg"
    },
    "Human Hemoglobin (complex)": {
        "sequence": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
        "image_url": "https://raw.githubusercontent.com/2-0-0-5-emil/protein-structure-app/main/human.jpeg"
    }
}

# ====== SIDEBAR CONTROLS ======
with st.sidebar:
    st.markdown("""
    <div style="text-align: center; margin-bottom: 0.7rem;">
        <h2 style="color: white; margin: 0;">Protein Input</h2>
    </div>
    """, unsafe_allow_html=True)
    
    example_sequences = {
        "Current sequence": "",
        "T4 Lysozyme (small)": protein_templates["T4 Lysozyme (small)"]["sequence"],
        "GFP (medium)": protein_templates["GFP (medium)"]["sequence"],
        "Human Hemoglobin (complex)": protein_templates["Human Hemoglobin (complex)"]["sequence"]
    }

    selected_example = st.selectbox(
        "Select an example protein", 
        options=list(example_sequences.keys()),
        key="protein_select"
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
        height=110,
        key="seq_input"
    )

    st.markdown("---", unsafe_allow_html=True)
    st.markdown('<h3 style="color:white;margin-bottom:0.4rem;">Visualization Options</h3>', unsafe_allow_html=True)
    
    color_scheme = st.selectbox(
        "Color Scheme", 
        ["spectrum", "chain", "residue", "secondary structure"],
        key="color_scheme"
    )
    
    spin = st.checkbox("Spin Structure", value=True, key="spin")

    predict_clicked = st.button("Predict Structure", key="predict_button", type="primary")
    reset = st.button("Reset", key="reset")

    if predict_clicked:
        st.session_state.prediction_made = True
        st.session_state.sequence = st.session_state.seq_input

    if reset:
        keys_to_clear = list(st.session_state.keys())
        for key in keys_to_clear:
            del st.session_state[key]
        st.session_state.prediction_made = False
        st.session_state.sequence = ""
        st.session_state.color_scheme = "spectrum"
        st.session_state.spin = True
        st.session_state.plddts = []

# ====== FUNCTION TO RUN PREDICTION AND SET STATE ======
def run_prediction_for_sequence(seq):
    st.session_state.sequence = seq
    st.session_state.prediction_made = True

# ====== HELPER FUNCTIONS ======
def render_mol(pdb, color_scheme='spectrum', spin=True):
    pdbview = py3Dmol.view()
    pdbview.addModel(pdb, 'pdb')
    pdbview.setStyle({'cartoon': {'color': color_scheme}})
    pdbview.setBackgroundColor('white')
    pdbview.zoomTo()
    pdbview.spin(spin)
    showmol(pdbview, height=500, width=800)

def get_confidence(plddt):
    if plddt >= 90:
        return "Very High"
    elif plddt >= 70:
        return "Confident"
    elif plddt >= 50:
        return "Low"
    else:
        return "Very Low"

def fetch_prediction(sequence):
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence)
    return response.content.decode('utf-8')

# ====== MAIN CONTENT AREA ======
if not st.session_state.prediction_made:
    st.markdown('<div class="main-content">', unsafe_allow_html=True)

    st.markdown("""
    <div class="custom-header">
        <h1 style="margin: 0; font-size:2.1rem; font-weight:800;">Protein Structure Predictor</h1>
        <p style="margin: 0.3rem 0 0; font-size: 1.05rem; opacity: 0.95;">
            Predict 3D protein structures using ESMFold's cutting-edge AI
        </p>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("### Try a preset protein sequence:")
    
    colA, colB, colC = st.columns(3)
    cols = [colA, colB, colC]
    for idx, (protein_name, protein_info) in enumerate(protein_templates.items()):
        with cols[idx]:
            st.markdown(f"""
            <div class="protein-card">
                <div class="protein-title">{protein_name}</div>
            </div>
            """, unsafe_allow_html=True)
            
            st.image(protein_info['image_url'], use_container_width=True)
            
            if st.button(f"Predict this protein", key=f"predict_{protein_name}_btn"):
                run_prediction_for_sequence(protein_info["sequence"])

    st.markdown('</div>', unsafe_allow_html=True)

else:
    sequence = st.session_state.sequence
    color_scheme = st.session_state.color_scheme
    spin = st.session_state.spin
    
    with st.spinner('Predicting protein structure...'):
        pdb_string = fetch_prediction(sequence)
        
        with open('predicted.pdb', 'w') as f:
            f.write(pdb_string)

        struct = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
        plddt = round(struct.b_factor.mean(), 2)
        seq_length = len(sequence)
        confidence = get_confidence(plddt)

        st.session_state.plddts = struct.b_factor.tolist()

        tab1, tab2, tab3 = st.tabs(["3D Structure", "Sequence Analysis", "PDB Data"])

        with tab1:
            st.subheader("Predicted Protein Structure")
            render_mol(pdb_string, color_scheme=color_scheme, spin=spin)

            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Sequence Length", seq_length)
            with col2:
                st.metric("Average pLDDT", plddt)
            with col3:
                st.metric("Confidence", confidence, 
                         help="pLDDT score indicates prediction confidence (0-100)")

        with tab2:
            st.subheader('Sequence Analysis')

            if st.session_state.plddts:
                st.markdown("""
                ### pLDDT per Residue
                pLDDT (predicted Local Distance Difference Test) estimates confidence per residue:
                - 90+ : Very high confidence
                - 70-90: High confidence
                - 50-70: Medium confidence
                - Below 50: Low confidence
                """)

                plddt_chart_data = pd.DataFrame({
                    "residue": list(range(1, len(st.session_state.plddts) + 1)), 
                    "pLDDT": st.session_state.plddts
                })
                st.line_chart(plddt_chart_data.set_index("residue"))

            if sequence:
                clean_seq = re.sub(r'[^A-Za-z]', '', sequence)
                aa_counts = {}
                for aa in clean_seq:
                    aa_counts[aa] = aa_counts.get(aa, 0) + 1

                st.markdown("### Amino Acid Composition")
                aa_data = pd.DataFrame({
                    "Amino Acid": list(aa_counts.keys()), 
                    "Count": list(aa_counts.values())
                })
                st.bar_chart(aa_data.set_index("Amino Acid"))

        with tab3:
            st.subheader("PDB Format Data")
            st.download_button(
                "Download PDB File", 
                pdb_string, 
                file_name="predicted.pdb", 
                mime="text/plain"
            )
            st.code(pdb_string, language='pdb')

# ====== FOOTER ======
st.markdown("""
<div class="custom-footer">
    <hr style="border: 0.5px solid #e2e8f0; margin: 1rem 0 0 0;">
    <p>Protein Structure Predictor v1.0 | Powered by ESMFold API</p>
    <p style="font-size: 0.8rem; margin-top: 0.25rem;">
        For research use only | Not for clinical or diagnostic use
    </p>
</div>
""", unsafe_allow_html=True)
