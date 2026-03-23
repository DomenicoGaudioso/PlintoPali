# -*- coding: utf-8 -*-
import json
import streamlit as st
from src import (
    DatiPlintoPali,
    DEFAULT_STRAT,
    valida_dati,
    calcola_plinto_pali,
    tabella_sintesi,
    tabella_pali,
    tabella_pali_comparativa,
    figura_geometria,
    figura_output,
    figura_comparativa,
    genera_warning,
    genera_note,
    export_json,
)

DEFAULTS = {
    'B': 4.0,
    'L': 4.0,
    'n_x': 2,
    'n_y': 2,
    'interasse_x': 2.5,
    'interasse_y': 2.5,
    'diametro_palo': 0.8,
    'lunghezza_palo': 18.0,
    'Nq': 25.0,
    'Nc': 9.0,
    'beta': 0.30,
    'alpha': 0.70,
    'gamma_sicurezza': 2.0,
    'N': 6000.0,
    'Mx': 800.0,
    'My': 600.0,
    'kh': 0.15,
    'kv': 0.05,
    'falda': 99.0,
    'stratigrafia_csv': DEFAULT_STRAT,
}

st.set_page_config(page_title='Plinto su pali', layout='wide')
st.title('PlintoPali - Statico, sismico e falda')
st.caption('v2.0: warning automatici, note tecniche, tabella comparativa, grafici avanzati.')

with st.sidebar:
    st.header('Import / Export input')
    up = st.file_uploader('Reimporta input JSON', type=['json'], key='plinto_json')
    defaults = DEFAULTS.copy()
    if up is not None:
        try:
            defaults.update(json.load(up))
            st.success('Input importati.')
        except Exception:
            st.error('JSON non valido.')

    st.header('Geometria plinto')
    B = st.number_input('Dimensione B [m]', 0.5, 50.0, float(defaults['B']), 0.1)
    L = st.number_input('Dimensione L [m]', 0.5, 50.0, float(defaults['L']), 0.1)
    n_x = st.number_input('Numero pali in x [-]', 1, 10, int(defaults['n_x']), 1)
    n_y = st.number_input('Numero pali in y [-]', 1, 10, int(defaults['n_y']), 1)
    interasse_x = st.number_input('Interasse pali in x [m]', 0.2, 10.0, float(defaults['interasse_x']), 0.1)
    interasse_y = st.number_input('Interasse pali in y [m]', 0.2, 10.0, float(defaults['interasse_y']), 0.1)
    diametro_palo = st.number_input('Diametro palo [m]', 0.1, 5.0, float(defaults['diametro_palo']), 0.05)
    lunghezza_palo = st.number_input('Lunghezza palo [m]', 1.0, 80.0, float(defaults['lunghezza_palo']), 0.5)

    st.header('Azioni')
    N = st.number_input('N [kN]', 0.0, 1e7, float(defaults['N']), 100.0)
    Mx = st.number_input('Mx [kNm]', -1e7, 1e7, float(defaults['Mx']), 50.0)
    My = st.number_input('My [kNm]', -1e7, 1e7, float(defaults['My']), 50.0)

    st.header('Capacità palo')
    Nq = st.number_input('Nq [-]', 1.0, 100.0, float(defaults['Nq']), 1.0)
    Nc = st.number_input('Nc [-]', 1.0, 20.0, float(defaults['Nc']), 0.5)
    beta = st.number_input('β attrito laterale [-]', 0.0, 2.0, float(defaults['beta']), 0.05)
    alpha = st.number_input('α aderenza [-]', 0.0, 2.0, float(defaults['alpha']), 0.05)
    gamma_sicurezza = st.number_input('γ sicurezza palo [-]', 1.1, 5.0, float(defaults['gamma_sicurezza']), 0.1)

    st.header('Sismica e falda')
    kh = st.number_input('kh [-]', 0.0, 1.0, float(defaults['kh']), 0.01)
    kv = st.number_input('kv [-]', 0.0, 1.0, float(defaults['kv']), 0.01)
    falda = st.number_input('Profondità falda [m]', 0.0, 100.0, float(defaults['falda']), 0.1)

    st.header('Stratigrafia')
    st.caption('Righe: spessore,gamma_dry,gamma_sat,phi,cu,k')
    stratigrafia_csv = st.text_area('Stratigrafia', value=str(defaults['stratigrafia_csv']), height=150)

# --- Costruzione dataclass e calcolo ---
d = DatiPlintoPali(
    B, L, int(n_x), int(n_y),
    interasse_x, interasse_y,
    diametro_palo, lunghezza_palo,
    Nq, Nc, beta, alpha,
    gamma_sicurezza, N, Mx, My,
    kh, kv, falda, stratigrafia_csv,
)
err = valida_dati(d)
if err:
    for e in err:
        st.error(e)
    st.stop()

r = calcola_plinto_pali(d)
df_sintesi = tabella_sintesi(r)
current = {k: v for k, v in d.__dict__.items()}
warnings = genera_warning(d, r)
note = genera_note(d, r)
n_tot = d.n_x * d.n_y

# --- Metriche ---
c1, c2, c3, c4, c5 = st.columns(5)
c1.metric('N° pali totale', f"{n_tot}")
c2.metric('Qamm palo [kN]', f"{r['Qamm_palo']:.0f}")
c3.metric('FS min statico', f"{float(min(r['statico']['FS'])):.2f}")
c4.metric('FS min sismico', f"{float(min(r['sismico']['FS'])):.2f}")

# Avvisi nelle metriche: icona se ci sono warning
if warnings:
    c5.metric('Avvisi', f"{len(warnings)} ⚠")
else:
    c5.metric('Avvisi', "Nessuno ✓")

# --- Tab ---
t1, t2, t3, t4, t5 = st.tabs([
    'Sintesi',
    'Stratigrafia',
    'Geometria Plotly',
    'Output Plotly',
    'Note tecniche',
])

with t1:
    st.subheader('Riepilogo globale')
    st.dataframe(df_sintesi, use_container_width=True)

    st.subheader('Tabella comparativa pali')
    df_comp = tabella_pali_comparativa(r)
    st.dataframe(df_comp, use_container_width=True)

    st.divider()
    col_a, col_b = st.columns(2)
    with col_a:
        st.download_button(
            'Salva input JSON',
            export_json(current),
            'plintopali_input.json',
            'application/json',
        )
    with col_b:
        st.download_button(
            'Scarica tabella comparativa CSV',
            df_comp.to_csv(index=False).encode('utf-8'),
            'plintopali_reazioni_comparativa.csv',
            'text/csv',
        )

with t2:
    st.subheader('Profilo stratigrafico')
    st.dataframe(r['stratigrafia'], use_container_width=True)

with t3:
    st.subheader('Pianta del plinto su pali')
    st.plotly_chart(figura_geometria(d, r), use_container_width=True)

with t4:
    st.subheader('Reazioni sui pali - caso statico')
    st.plotly_chart(figura_output(r, 'statico'), use_container_width=True)
    st.subheader('Reazioni sui pali - caso sismico')
    st.plotly_chart(figura_output(r, 'sismico'), use_container_width=True)
    st.subheader('Confronto statico vs sismico')
    st.plotly_chart(figura_comparativa(r), use_container_width=True)

with t5:
    st.subheader('Avvertimenti di progetto')
    if warnings:
        for w in warnings:
            if w.startswith('⚠'):
                st.warning(w)
            else:
                st.info(w)
    else:
        st.success('Nessun avvertimento: il progetto soddisfa tutti i controlli automatici.')

    st.subheader('Note tecniche')
    for nota in note:
        st.markdown(f"- {nota}")
