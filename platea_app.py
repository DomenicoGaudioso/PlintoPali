# -*- coding: utf-8 -*-
import json
import pandas as pd
import streamlit as st
from platea_src import (
    DatiPlatea,
    valida_dati_platea,
    calcola_platea_fem,
    figura_geometria_platea,
    figura_risultati_platea,
)

DEFAULTS = {
    'B': 10.0, 'L': 12.0, 'spessore': 0.8, 'E_cls_MPa': 30000.0,
    'k_winkler_kPa_m': 15000.0, 'mesh_size': 1.0,
}

DEFAULT_COLUMNS = pd.DataFrame([
    {"x": 2.0, "y": 2.0, "P_kN": 1200.0, "Mx_kNm": 50.0, "My_kNm": -30.0},
    {"x": 8.0, "y": 2.0, "P_kN": 1500.0, "Mx_kNm": 0.0, "My_kNm": 0.0},
    {"x": 2.0, "y": 10.0, "P_kN": 1100.0, "Mx_kNm": 0.0, "My_kNm": 0.0},
    {"x": 8.0, "y": 10.0, "P_kN": 1600.0, "Mx_kNm": -40.0, "My_kNm": 25.0},
])

st.set_page_config(page_title='PlateaFEM', layout='wide')
st.title('PlateaFEM - Analisi di Platee di Fondazione')
st.markdown("Strumento per l'analisi di platee di fondazione su suolo alla Winkler con metodo degli Elementi Finiti (FEM) basato su **OpenSeesPy**.")

with st.sidebar:
    st.header('📐 Geometria Platea')
    B = st.number_input('Dimensione B (dir. X) [m]', 1.0, 100.0, float(DEFAULTS['B']), 0.5)
    L = st.number_input('Dimensione L (dir. Y) [m]', 1.0, 100.0, float(DEFAULTS['L']), 0.5)
    spessore = st.number_input('Spessore Platea H [m]', 0.2, 5.0, float(DEFAULTS['spessore']), 0.1)
    
    st.header('⚙️ Proprietà Materiali e Terreno')
    E_cls_MPa = st.number_input('Modulo E Cls [MPa]', 10000.0, 60000.0, float(DEFAULTS['E_cls_MPa']), 1000.0)
    k_winkler_kPa_m = st.number_input('Modulo di Winkler ks [kPa/m]', 1000.0, 500000.0, float(DEFAULTS['k_winkler_kPa_m']), 1000.0)

    st.header('FEM')
    mesh_size = st.number_input('Dimensione Mesh [m]', 0.2, 5.0, float(DEFAULTS['mesh_size']), 0.1)

    st.header('🏛️ Carichi Pilastri')
    st.caption('Definire posizione e carichi dei pilastri. Le coordinate (0,0) sono nell\'angolo in basso a sinistra.')
    
    edited_columns_df = st.data_editor(
        DEFAULT_COLUMNS,
        num_rows="dynamic",
        use_container_width=True,
        column_config={
            "x": st.column_config.NumberColumn("x [m]", format="%.2f"),
            "y": st.column_config.NumberColumn("y [m]", format="%.2f"),
            "P_kN": st.column_config.NumberColumn("P [kN]", format="%.1f"),
            "Mx_kNm": st.column_config.NumberColumn("Mx [kNm]", format="%.1f"),
            "My_kNm": st.column_config.NumberColumn("My [kNm]", format="%.1f"),
        },
        key="columns_editor"
    )

# Costruzione istanza Dati
dati_platea = DatiPlatea(
    B=B, L=L, spessore=spessore, E_cls_MPa=E_cls_MPa,
    k_winkler_kPa_m=k_winkler_kPa_m, mesh_size=mesh_size,
    pilastri_df=edited_columns_df
)

# Validazione
err = valida_dati_platea(dati_platea)
if err:
    for e in err:
        st.error(e)
    st.stop()

# Calcolo
try:
    with st.spinner("Analisi FEM della platea in corso..."):
        risultati = calcola_platea_fem(dati_platea)

    cedimento_max = risultati['cedimenti_mm'].max().max()
    cedimento_min = risultati['cedimenti_mm'].min().min()
    pressione_max = risultati['pressioni_kPa'].max().max()

    c1, c2, c3 = st.columns(3)
    c1.metric("Cedimento Max [mm]", f"{cedimento_max:.2f}")
    c2.metric("Cedimento Min [mm]", f"{cedimento_min:.2f}")
    c3.metric("Pressione Max [kPa]", f"{pressione_max:.1f}")

    # Tabs per i risultati
    t1, t2, t3, t4 = st.tabs([
        '📐 Geometria e Mesh', 
        '📉 Cedimenti', 
        '📈 Momenti Flettenti',
        '📊 Pressioni sul Terreno'
    ])

    with t1:
        st.subheader("Geometria, Mesh FEM e Posizione Carichi")
        st.plotly_chart(figura_geometria_platea(dati_platea, risultati), use_container_width=True)

    with t2:
        st.subheader("Mappa dei Cedimenti Verticali [mm]")
        st.plotly_chart(figura_risultati_platea(dati_platea, risultati, 'cedimenti_mm', 'Cedimenti [mm]'), use_container_width=True)

    with t3:
        st.subheader("Mappe dei Momenti Flettenti [kNm/m]")
        col1, col2 = st.columns(2)
        with col1:
            st.plotly_chart(figura_risultati_platea(dati_platea, risultati, 'Mxx_kNm_m', 'Momento Mxx [kNm/m]'), use_container_width=True)
        with col2:
            st.plotly_chart(figura_risultati_platea(dati_platea, risultati, 'Myy_kNm_m', 'Momento Myy [kNm/m]'), use_container_width=True)
        st.plotly_chart(figura_risultati_platea(dati_platea, risultati, 'Mxy_kNm_m', 'Momento Mxy [kNm/m]'), use_container_width=True)

    with t4:
        st.subheader("Mappa delle Pressioni sul Terreno [kPa]")
        st.plotly_chart(figura_risultati_platea(dati_platea, risultati, 'pressioni_kPa', 'Pressione sul terreno [kPa]'), use_container_width=True)

except Exception as e:
    st.error(f"Errore critico durante l'analisi FEM: {e}")
    st.exception(e)

```
```diff