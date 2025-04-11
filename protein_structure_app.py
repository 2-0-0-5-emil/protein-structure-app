import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import re
import io
import biotite.structure.io as bsio

# Streamlit Page Setup
st.set_page_config(layout="wide", page_title="Protein Structure Predictor")

st.sidebar.title('🧬 ESMFold Protein Predictor')
st.sidebar.write('Enter your protein sequence below and click **Predict** to view its 3D structure and pLDDT score.')

DEFAULT_SEQ = "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ"
sequence_input = st.sidebar.text_area('Input Protein Sequence', DEFAULT_SEQ, height=250)

color_option = st.sidebar.selectbox("Visualization Color Scheme", ["spectrum", "chain", "sstruc", "residue"])

# Function: Validate input
def is_valid_protein_sequence(seq):
    clean_seq = re.sub(r'[^A-Za-z]', '', seq)
    valid_aa = set("ACDEFGHIKLMNPQRSTVWYBJOUXZ")
    return all(aa.upper() in valid_aa for aa in clean_seq) and len(clean_seq) > 0

# Function: Render 3D mol
def render_mol(pdb_data, color_scheme="spectrum"):
    view = py3Dmol.view(width=800, height=500)
    view.addModel(pdb_data, 'pdb')
    if color_scheme == "spectrum":
        view.setStyle({'cartoon': {'color': 'spectrum'}})
    elif color_scheme == "chain":
        view.setStyle({'cartoon': {'colorscheme': 'chain'}})
    elif color_scheme == "sstruc":
        view.setStyle({'cartoon': {'colorscheme': 'sstruc'}})
    elif color_scheme == "residue":
        view.setStyle({'cartoon': {'colorscheme': 'amino'}})
    else:
        view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.setBackgroundColor('white')
    view.zoomTo()
    view.spin(True)
    showmol(view, height=500, width=800)

# Function: Extract plDDT and metrics
def extract_metrics(pdb_string):
    pdb_file = io.StringIO(pdb_string)
    struct = bsio.load_structure(pdb_file, extra_fields=["b_factor"])
    b_value = round(struct.b_factor.mean(), 2)
    length = len(set(struct.res_id))
    return b_value, length

# Main predict function
def predict_structure(seq):
    with st.spinner("Predicting structure using ESMFold..."):
        headers = {'Content-Type': 'application/x-www-form-urlencoded'}
        response = requests.post("https://api.esmatlas.com/foldSequence/v1/pdb/", headers=headers, data=seq)
        pdb_data = response.content.decode('utf-8')
        b_score, length = extract_metrics(pdb_data)
        
        st.success("✅ Prediction complete!")
        st.subheader("🧬 3D Visualization of Protein Structure")
        render_mol(pdb_data, color_scheme=color_option)

        st.subheader("📈 Metrics")
        st.markdown(f"- **Protein Length**: {length} residues  
                     - **Average plDDT Score**: `{b_score}` (Confidence out of 100)")

        st.download_button("📥 Download PDB", pdb_data, file_name="predicted.pdb")

# Predict Button
if st.sidebar.button("🔍 Predict Structure"):
    if is_valid_protein_sequence(sequence_input):
        predict_structure(sequence_input)
    else:
        st.error("⚠️ Please enter a valid protein sequence using only amino acid codes (A-Z).")

else:
    st.warning("👈 Paste a valid protein sequence and hit Predict!")
