import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio

# --- Config ---
st.set_page_config(layout='wide')

# --- Helper Functions ---
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

# --- UI ---
st.title('🔬 Protein Structure Predictor using ESMFold')
st.markdown('---')

# Sidebar Inputs
st.sidebar.title('Protein Input & Settings')
def_seq = "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ"

example_btn = st.sidebar.button("Use Example Protein")
if example_btn:
    st.session_state.sequence = def_seq
sequence = st.sidebar.text_area("Enter Protein Sequence", st.session_state.get("sequence", ""), height=250)

# Visualization settings
st.sidebar.title("🔧 Visualization Options")
color_scheme = st.sidebar.selectbox("Color Scheme", ["spectrum", "chain", "residue", "secondary structure"])
spin = st.sidebar.checkbox("Spin Structure", value=True)
reset = st.sidebar.button("Reset", type="primary")
if reset:
    st.session_state.clear()
    st.experimental_rerun()

predict = st.sidebar.button("Predict Structure")

# Main App Tabs
if predict and sequence:
    pdb_string = fetch_prediction(sequence)
    with open('predicted.pdb', 'w') as f:
        f.write(pdb_string)

    struct = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
    plddt = round(struct.b_factor.mean(), 2)
    seq_length = len(sequence)
    confidence = get_confidence(plddt)

    tab1, tab2, tab3 = st.tabs(["3D Structure", "Sequence Analysis", "PDB Data"])

    with tab1:
        st.subheader("Predicted Protein Structure")
        render_mol(pdb_string, color_scheme=color_scheme, spin=spin)
        
        col1, col2, col3 = st.columns(3)
        col1.metric("Sequence Length", seq_length)
        col2.metric("Average pLDDT", plddt)
        col3.metric("Confidence", confidence)

    with tab2:
        st.subheader("Sequence Analysis")
        st.code(sequence)

    with tab3:
        st.subheader("PDB Format Data")
        st.download_button("Download PDB", pdb_string, file_name="predicted.pdb", mime="text/plain")
        st.text_area("PDB Data", pdb_string, height=300)

elif not predict:
    st.info("👈 Paste a sequence and click Predict to begin.")
