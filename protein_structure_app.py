import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio
import re
import pandas as pd
import random
st.set_page_config(layout='wide')
st.markdown(
    """
    <style>
    .stApp {
        background-image: url("https://imgs.search.brave.com/-R07OhwFpq0uHMugmpZQylZSCQQuc8wjus46YHg8640/rs:fit:500:0:0:0/g:ce/aHR0cHM6Ly90My5m/dGNkbi5uZXQvanBn/LzAyLzgzLzQ2Lzcw/LzM2MF9GXzI4MzQ2/NzA1NV8xeWxVeE1q/dkE3bG9PWW5rUlZ0/bnlqM2NxajZuWEZW/UC5qcGc");
        background-size: cover;
        background-position: center;
        background-attachment: fixed;
    }
    [data-testid="stSidebar"] > div:first-child {
        background-image: url("https://i.pinimg.com/474x/74/ee/c2/74eec2a08cc9dd5328d78a869e965884.jpg");
        background-size: cover;
        background-position: center;
        background-repeat: no-repeat;
    }
    [data-testid="stSidebar"] {
        color: white;
    }
    section[data-testid="stSidebar"] h1, section[data-testid="stSidebar"] h2, section[data-testid="stSidebar"] h3 {
        color: white;
    }
    </style>
    """,
    unsafe_allow_html=True
)
st.sidebar.title("Navigation")
page = st.sidebar.selectbox("Choose a page:", ["Protein Quiz", "Protein Structure Predictor"])
if page == "Protein Quiz":
    st.markdown('<div class="quiz-box">', unsafe_allow_html=True)
    st.title("🧪 Mini Protein Quiz")
    st.markdown("Test your knowledge about proteins!")
    st.markdown("### How many essential amino acids must be obtained from the diet?")
    answer1 = st.radio("Choose one:", ["5", "10", "20", "22"], key="q1")
    st.markdown("### What are the building blocks of proteins?")
    answer2 = st.radio("Choose one:", ["Nucleotides", "Fatty acids", "Amino acids", "Monosaccharides"], key="q2")
    st.markdown("### How many standard amino acids are there?")
    answer3 = st.radio("Choose one:", ["10", "12", "20", "100"], key="q3")
    st.markdown("### Which organelle in the cell helps make proteins?")
    answer4 = st.radio("Choose one:", ["Mitochondria", "Ribosome", "Nucleus", "Golgi Apparatus"], key="q4")
    if st.button("Submit Answers"):
        score = 0
        if answer1 == "10":
            score += 1
        if answer2 == "Amino acids":
            score += 1
        if answer3 == "20":
            score += 1
        if answer4 == "Ribosome":
            score += 1
        st.markdown(f"### 🎯 Your Score: {score}/4")

        if score == 4:
            st.success("✅ Perfect! You’re a protein pro! 🎉")
            st.balloons()
            st.markdown("""
            <audio autoplay>
                <source src="https://www.soundjay.com/buttons/sounds/button-4.mp3" type="audio/mpeg">
            </audio>
            """, unsafe_allow_html=True)
        else:
            st.error("❌ Some answers are incorrect. Try again!")
            st.markdown("""
            <audio autoplay>
                <source src="https://www.soundjay.com/button/beep-07.wav" type="audio/wav">
            </audio>
            """, unsafe_allow_html=True)

    st.markdown('</div>', unsafe_allow_html=True)

# --- PAGE 2: Protein Structure Predictor ---
elif page == "Protein Structure Predictor":
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

    st.title('🔬 Protein Structure Predictor using ESMFold')

    st.sidebar.title('Protein Input & Settings')

    col_left, col_right = st.columns([1, 3])

    with col_left:
        st.image("https://i.pinimg.com/474x/79/7a/c5/797ac574deece2d7b8e2ef844a4fe571.jpg", width=200)

    with col_right:
        if st.button("Click the Mascot for a Fun Protein Fact!"):
            st.success("""
            **Fun Protein Fact**: Proteins are essential for almost every function in our bodies!
            They are responsible for building tissues, fighting infections, and much more.

            **Did you know?** A protein can have over 1,000 amino acids in its structure!
            """)
        else:
            st.markdown("_Click the mascot to learn something cool!_")

    st.markdown("---")

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

    st.sidebar.title("🔧 Visualization Options")
    color_scheme = st.sidebar.selectbox("Color Scheme", ["spectrum", "chain", "residue", "secondary structure"])
    spin = st.sidebar.checkbox("Spin Structure", value=True)
    reset = st.sidebar.button("Reset", type="primary")
    if reset:
        for key in st.session_state.keys():
            del st.session_state[key]
        st.rerun()

    predict = st.sidebar.button("Predict Structure")

    if predict and sequence:
        gif_placeholder = st.empty()
        gif_placeholder.markdown(
        """
        <div style="display: flex; justify-content: center; margin: 20px 0;">
            <img src="https://media.tenor.com/vp3V50Hs-B8AAAAm/loading-waiting.webp">
        </div>
        """,
        unsafe_allow_html=True
    )

        with st.spinner("Predicting Structure..."):
            pdb_string = fetch_prediction(sequence)

        gif_placeholder.empty()
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

            if st.session_state.plddts:
                st.write("### pLDDT per Residue")
                plddt_chart_data = pd.DataFrame({"residue": list(range(1, len(st.session_state.plddts) + 1)), "pLDDT": st.session_state.plddts})
                st.line_chart(plddt_chart_data.set_index("residue"))

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
