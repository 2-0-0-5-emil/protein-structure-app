import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio
import matplotlib.pyplot as plt
import re

# UI Setup
st.sidebar.title('🔬 Protein Structure Predictor')
st.sidebar.write("Enter a protein sequence to visualize its predicted 3D structure using ESMFold.")

# Protein sequence input
DEFAULT_SEQ = "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ"
txt = st.sidebar.text_area('Input Sequence', DEFAULT_SEQ, height=275)

# Color scheme selector
color_scheme = st.sidebar.selectbox("Color Scheme", ["spectrum", "chain", "sstruc", "residue"])

# Validate sequence
def is_valid_sequence(seq):
    return re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWYBJOUXZ]+', seq.strip().upper()) is not None

# Visualize molecule
def render_mol(pdb, scheme="spectrum"):
    view = py3Dmol.view(width=800, height=500)
    view.addModel(pdb, 'pdb')
    
    color_map = {
        "spectrum": {"cartoon": {"color": "spectrum"}},
        "chain": {"cartoon": {"color": "chain"}},
        "sstruc": {"cartoon": {"color": "sstruc"}},
        "residue": {"cartoon": {"color": "amino"}},
    }
    
    view.setStyle(color_map.get(scheme, {"cartoon": {"color": "spectrum"}}))
    view.setBackgroundColor('white')
    view.zoomTo()
    view.zoom(2, 800)
    view.spin(True)
    showmol(view, height=500, width=800)

# Predict and visualize
def update(sequence=txt):
    if not is_valid_sequence(sequence):
        st.error("❌ Invalid protein sequence! Use only valid amino acid codes (A-Z, no numbers or symbols).")
        return

    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence)
    pdb_string = response.content.decode('utf-8')

    with open('predicted.pdb', 'w') as f:
        f.write(pdb_string)

    struct = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
    b_value = round(struct.b_factor.mean(), 4)

    st.subheader('🧬 3D Structure')
    render_mol(pdb_string, scheme=color_scheme)

    st.subheader('📊 plDDT Score')
    st.info(f'Average plDDT: {b_value}')
    st.caption("plDDT (predicted Local Distance Difference Test) scores indicate prediction confidence (0–100).")

    st.subheader("📈 Per-residue plDDT Scores")
    plddts = struct.b_factor.tolist()
    fig, ax = plt.subplots()
    ax.plot(plddts, color='green')
    ax.set_xlabel("Residue Index")
    ax.set_ylabel("plDDT Score")
    st.pyplot(fig)

    st.download_button(
        label="Download PDB File",
        data=pdb_string,
        file_name='predicted.pdb',
        mime='text/plain',
    )

# Prediction trigger
if st.sidebar.button('Predict'):
    update(txt)
else:
    st.warning('👈 Enter a valid protein sequence and click Predict.')
