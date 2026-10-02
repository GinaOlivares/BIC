"""Application entrypoint with explicit public navigation labels."""
import streamlit as st

st.set_page_config(page_title='BIC Explorer', layout='wide')
page = st.navigation([
    st.Page('pages/0_BIC_Explorer.py', title='BIC Explorer', default=True),
    st.Page('pages/1_Waveguide_Modes.py', title='Waveguide Modes', url_path='Waveguide_Modes'),
])
page.run()
