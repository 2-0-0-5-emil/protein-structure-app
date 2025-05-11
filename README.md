🔬 Protein Structure Prediction App
Live Demo: Try it on Streamlit! 🚀
![Screenshot 2025-05-11 210306](https://github.com/user-attachments/assets/34dc9e63-a883-4b1b-a4e5-1b1cc7179826)

📌 About the Project
This Streamlit web application allows users to:

🧬 Predict protein structures from amino acid sequences.

🔍 Visualize 3D structures interactively using Py3Dmol.

🧪 Explore examples and adjust visualization settings.

Built to explore the intersection of AI and bioinformatics, this project uses AlphaFold2 (via BioColabFold) under the hood.

⚙️ Features
🧠 Structure prediction using a simplified AlphaFold2 pipeline.

🖱️ Interactive 3D viewer with color, spin, and reset controls.

📂 Upload or paste a sequence in FASTA format.

🧪 Try out example proteins in one click.

🚀 Getting Started

git clone https://github.com/2-0-0-5-emil/protein-structure-app.git
cd protein-structure-app
pip install -r requirements.txt
streamlit run app.py

📁 File Structure

protein-structure-app/
│
├── app.py                 # Main Streamlit app
├── utils.py               # Utility functions for prediction
├── assets/                # Icons and screenshots
├── requirements.txt       # Python dependencies
└── README.md              # Project documentation


🛠️ Built With
Streamlit

ColabFold

Py3Dmol

BioPython

🙋‍♂️ Author
Made with 🧠 and ❤️ by Emil Jinu
LinkedIn | GitHub

⭐ Contribute
If you like this project, consider giving it a ⭐ and contributing! Open issues or suggestions are always welcome.

