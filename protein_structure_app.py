import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio

st.sidebar.title('🔬 Protein Structure Predictor')
st.sidebar.markdown(
    "Powered by [ESMFold](https://esmatlas.com/about), this app predicts and visualizes protein structures directly from amino acid sequences.\n\n"
    "Enter a valid protein sequence and click **Predict** to generate a 3D model and confidence score (plDDT).\n\n"
    "Built with ❤️ using Streamlit, py3Dmol, and Biotite."
)

def render_mol(pdb):
    pdbview = py3Dmol.view()
    pdbview.addModel(pdb, 'pdb')
    pdbview.setStyle({'cartoon': {'color': 'spectrum'}})
    pdbview.setBackgroundColor('white')
    pdbview.zoomTo()
    pdbview.zoom(2, 800)
    pdbview.spin(True)
    showmol(pdbview, height=500, width=800)

DEFAULT_SEQ = "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ"
txt = st.sidebar.text_area('Input sequence', DEFAULT_SEQ, height=275)

def update(sequence=txt):
    headers = {
        'Content-Type': 'application/x-www-form-urlencoded',
    }
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence)
    name = sequence[:3] + sequence[-3:]
    pdb_string = response.content.decode('utf-8')

    with open('predicted.pdb', 'w') as f:
        f.write(pdb_string)

    struct = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
    b_value = round(struct.b_factor.mean(), 4)

    st.subheader('🧬 Visualization of Predicted Protein Structure')
    render_mol(pdb_string)

    st.subheader('📊 Prediction Confidence (plDDT)')
    st.write("plDDT is a per-residue confidence score for predicted structure accuracy (scale: 0–100).")
    st.info(f'Average plDDT: {b_value}')

    st.download_button(
        label="📥 Download PDB File",
        data=pdb_string,
        file_name='predicted.pdb',
        mime='text/plain',
    )

predict = st.sidebar.button('Predict', on_click=update)

if not predict:
    st.warning('👈 Enter a protein sequence and click Predict!')
