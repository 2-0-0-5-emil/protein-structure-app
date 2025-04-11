import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio
import re
import pandas as pd

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

example_sequences = {
    "Current sequence": "",
    "T4 Lysozyme (small)": "MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL",
    "GFP (medium)": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
    "Human Hemoglobin (complex)": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
}

selected_example = st.sidebar.selectbox("Select an example protein", options=list(example_sequences.keys()))
sequence_input = example_sequences[selected_example]

if selected_example == "Current sequence":
    sequence_input = st.session_state.get("sequence", "")
    st.session_state.sequence = sequence_input
else:
    st.session_state.sequence = sequence_input

sequence = st.sidebar.text_area("Enter Protein Sequence", sequence_input, height=250)

# Visualization settings
st.sidebar.title("🔧 Visualization Options")
color_scheme = st.sidebar.selectbox("Color Scheme", ["spectrum", "chain", "residue", "secondary structure"])
spin = st.sidebar.checkbox("Spin Structure", value=True)
reset = st.sidebar.button("Reset", type="primary")
if reset:
    for key in st.session_state.keys():
        del st.session_state[key]
    st.rerun()

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

    st.session_state.plddts = struct.b_factor.tolist()

    tab1, tab2, tab3 = st.tabs(["3D Structure", "Sequence Analysis", "PDB Data"])

    with tab1:
        st.subheader("Predicted Protein Structure")
        render_mol(pdb_string, color_scheme=color_scheme, spin=spin)

        col1, col2, col3 = st.columns(3)
        col1.metric("Sequence Length", seq_length)
        col2.metric("Average pLDDT", plddt)
        col3.metric("Confidence", confidence)

    with tab2:
        st.subheader('Sequence Analysis')

        # Display pLDDT chart
        if st.session_state.plddts:
            st.write("### pLDDT per Residue")
            st.write("pLDDT (predicted Local Distance Difference Test) is a per-residue estimate of the confidence in prediction on a scale from 0-100.")
            st.write("- 90+ : Very high confidence")
            st.write("- 70-90: High confidence")
            st.write("- 50-70: Medium confidence")
            st.write("- Below 50: Low confidence")

            # Plot pLDDT values
            plddt_chart_data = pd.DataFrame({"residue": list(range(1, len(st.session_state.plddts) + 1)), "pLDDT": st.session_state.plddts})
            st.line_chart(plddt_chart_data.set_index("residue"))

        # Display amino acid composition
        if sequence:
            clean_seq = re.sub(r'[^A-Za-z]', '', sequence)
            aa_counts = {}
            for aa in clean_seq:
                aa_counts[aa] = aa_counts.get(aa, 0) + 1

            st.write("### Amino Acid Composition")
            aa_data = pd.DataFrame({"Amino Acid": list(aa_counts.keys()), "Count": list(aa_counts.values())})
            st.bar_chart(aa_data.set_index("Amino Acid"))

    with tab3:
        st.subheader("PDB Format Data")
        st.download_button("Download PDB", pdb_string, file_name="predicted.pdb", mime="text/plain")
        st.text_area("PDB Data", pdb_string, height=300)

elif not predict:
    st.info("👈 Paste a sequence and click Predict to begin.")
