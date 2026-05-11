import json
import streamlit as st

from src import (
    DEFAULT_STRAT,
    DatiPalo,
    valida_dati,
    calcola_palo,
    export_json,
    tabella_sintesi,
    tabella_domanda_capacita,
    valuta_esiti,
    tabella_confronto_palo,
    tabella_tensioni_effettive,
    figura_geometria,
    figura_spostamenti,
    figura_confronto_palo,
    genera_warning,
    genera_note,
    riepilogo_progetto,
)

# Using the same CSS from Muro app for consistency
PREMIUM_CSS = """
<style>
  :root {
    --ink: #13202b;
    --muted: #5b6774;
    --line: rgba(19, 32, 43, 0.12);
    --paper: #f6f7f8;
    --panel: #ffffff;
    --accent: #155e75;
    --ok: #166534;
    --warn: #a16207;
    --bad: #b91c1c;
  }
  .stApp {
    background: linear-gradient(180deg, #f7f8fa 0%, #eef2f5 100%);
  }
  .block-container {
    padding-top: 1.4rem;
    padding-bottom: 2rem;
  }
  .hero-shell {
    border: 1px solid rgba(255,255,255,.22);
    background: linear-gradient(135deg, #13202b 0%, #164e63 100%);
    color: white;
    padding: 24px 28px;
    border-radius: 14px;
    box-shadow: 0 18px 42px rgba(19, 32, 43, .16);
    margin-bottom: 1rem;
  }
  .hero-shell h1 {
    margin: 0 0 .35rem;
    font-size: 2.05rem;
    line-height: 1.05;
    letter-spacing: 0;
  }
  .hero-shell p {
    margin: 0;
    max-width: 860px;
    color: rgba(255,255,255,.82);
    font-size: 1.02rem;
  }
  .hero-badges {
    display: flex;
    gap: .5rem;
    flex-wrap: wrap;
    margin-top: 1rem;
  }
  .hero-badge {
    border: 1px solid rgba(255,255,255,.18);
    background: rgba(255,255,255,.10);
    border-radius: 8px;
    padding: .28rem .72rem;
    font-size: .78rem;
  }
  .section-card {
    background: var(--panel);
    border: 1px solid var(--line);
    border-radius: 10px;
    padding: 1rem 1.15rem;
    box-shadow: 0 8px 24px rgba(20, 39, 55, .05);
    margin-bottom: .85rem;
    color: var(--ink);
  }
  .section-card:empty {
    display: none;
  }
  .section-card h1,
  .section-card h2,
  .section-card h3,
  .section-card p,
  .section-card li,
  .section-card span {
    color: var(--ink);
  }
  .mini-grid {
    display: grid;
    grid-template-columns: repeat(3, minmax(0, 1fr));
    gap: .85rem;
    margin: 1rem 0 0;
  }
  .mini-card {
    background: rgba(255,255,255,.12);
    border: 1px solid rgba(255,255,255,.14);
    border-radius: 10px;
    padding: .9rem 1rem;
  }
  .mini-card .label {
    font-size: .75rem;
    text-transform: uppercase;
    letter-spacing: .08em;
    color: rgba(255,255,255,.68);
    margin-bottom: .3rem;
  }
  .mini-card .value {
    font-size: 1.28rem;
    font-weight: 700;
  }
  div[data-testid="stMetric"] {
    background: #ffffff;
    border: 1px solid var(--line);
    border-radius: 10px;
    padding: .8rem 1rem;
    box-shadow: 0 8px 20px rgba(20, 39, 55, .04);
  }
  div[data-testid="stMetric"] * {
    color: var(--ink) !important;
  }
  div[data-testid="stMetricLabel"] {
    color: var(--muted) !important;
  }
  div[data-testid="stMetricValue"] {
    color: var(--ink) !important;
  }
  .stTabs [data-baseweb="tab-list"] {
    gap: .35rem;
  }
  .stTabs [data-baseweb="tab"] {
    background: rgba(255,255,255,.62);
    border: 1px solid var(--line);
    border-radius: 8px;
    height: 42px;
    padding: 0 1rem;
  }
  .stTabs [data-baseweb="tab"] p {
    color: var(--ink);
  }
  .stTabs [aria-selected="true"] {
    background: #155e75;
    color: white;
    border-color: transparent;
  }
  .stTabs [aria-selected="true"] p {
    color: white;
  }
  .status-pill {
    display: inline-flex;
    align-items: center;
    border-radius: 999px;
    padding: .16rem .55rem;
    font-size: .76rem;
    font-weight: 700;
  }
  .status-ok { background: #dcfce7; color: var(--ok); }
  .status-warn { background: #fef3c7; color: var(--warn); }
  .status-ko { background: #fee2e2; color: var(--bad); }
  [data-testid="stMarkdownContainer"] h2,
  [data-testid="stMarkdownContainer"] h3,
  [data-testid="stMarkdownContainer"] p,
  [data-testid="stMarkdownContainer"] li {
    color: var(--ink);
  }
  [data-testid="stCaptionContainer"],
  [data-testid="stCaptionContainer"] * {
    color: var(--muted) !important;
  }
  .hero-shell h1,
  .hero-shell p,
  .hero-badge,
  .hero-badge * {
    color: white !important;
  }
</style>
"""

DEFAULTS = {
    "diametro": 0.6,
    "lunghezza": 15.0,
    "Nq": 20.0,
    "Nc": 9.0,
    "beta": 0.3,
    "alpha": 0.5,
    "N": 800.0,
    "H": 50.0,
    "E": 3.0e7,
    "I": 0.00636,
    "kh": 0.1,
    "kv": 0.05,
    "falda": 5.0,
    "n_elementi": 40,
    "stratigrafia_csv": DEFAULT_STRAT,
}

st.set_page_config(page_title="Palo Singolo", layout="wide")
st.markdown(PREMIUM_CSS, unsafe_allow_html=True)

st.markdown(
    """
    <section class="hero-shell">
      <h1>Palo Singolo</h1>
      <p>Analisi geotecnica di un palo singolo soggetto a carichi assiali e laterali, in condizioni statiche e pseudo-statiche.</p>
      <div class="hero-badges">
        <span class="hero-badge">Capacità Assiale</span>
        <span class="hero-badge">Risposta Laterale</span>
        <span class="hero-badge">Winkler (FD/FEM)</span>
        <span class="hero-badge">Curve p-y</span>
      </div>
      <div class="mini-grid">
        <div class="mini-card">
          <div class="label">Verifiche</div>
          <div class="value">SLU + SLE</div>
        </div>
        <div class="mini-card">
          <div class="label">Motore FEM</div>
          <div class="value">Interno</div>
        </div>
        <div class="mini-card">
          <div class="label">Normative</div>
          <div class="value">NTC2018 / EC7</div>
        </div>
      </div>
    </section>
    """,
    unsafe_allow_html=True,
)

with st.sidebar:
    st.header("Import / Export")
    uploaded_input = st.file_uploader("Importa input JSON", type=["json"], key="palo_json")
    defaults = DEFAULTS.copy()
    if uploaded_input is not None:
        try:
            defaults.update(json.load(uploaded_input))
            st.success("Input importati correttamente.")
        except Exception:
            st.error("JSON non valido.")

    st.header("Geometria e Materiale Palo")
    c1, c2 = st.columns(2)
    with c1:
        diametro = st.number_input("Diametro [m]", 0.2, 3.0, float(defaults["diametro"]), 0.1)
        E = st.number_input("Modulo E [kPa]", 1e6, 5e7, float(defaults["E"]), 1e6)
    with c2:
        lunghezza = st.number_input("Lunghezza [m]", 1.0, 50.0, float(defaults["lunghezza"]), 0.5)
        I = st.number_input("Inerzia I [m4]", 1e-4, 1.0, float(defaults["I"]), 1e-4, format="%.5f")

    st.header("Carichi in Testa")
    c1, c2 = st.columns(2)
    with c1:
        N = st.number_input("Carico Assiale N [kN]", 0.0, 10000.0, float(defaults["N"]), 50.0)
    with c2:
        H = st.number_input("Carico Laterale H [kN]", 0.0, 1000.0, float(defaults["H"]), 10.0)

    st.header("Parametri Geotecnici Assiali")
    c1, c2, c3 = st.columns(3)
    with c1:
        Nq = st.number_input("Nq [-]", 1.0, 100.0, float(defaults["Nq"]), 1.0)
        beta = st.number_input("beta [-]", 0.1, 1.0, float(defaults["beta"]), 0.05)
    with c2:
        Nc = st.number_input("Nc [-]", 5.0, 20.0, float(defaults["Nc"]), 0.5)
        alpha = st.number_input("alpha [-]", 0.1, 1.0, float(defaults["alpha"]), 0.05)
    with c3:
        falda = st.number_input("Profondità falda [m]", 0.0, 100.0, float(defaults["falda"]), 0.5)

    st.header("Sismica")
    c1, c2 = st.columns(2)
    with c1:
        kh = st.number_input("kh [-]", 0.0, 1.0, float(defaults["kh"]), 0.01)
    with c2:
        kv = st.number_input("kv [-]", 0.0, 1.0, float(defaults["kv"]), 0.01)

    st.header("Modello Numerico")
    n_elementi = st.number_input("Numero elementi (FD/FEM)", 10, 100, int(defaults["n_elementi"]), 5)
    usa_py = st.checkbox("Abilita curve p-y non lineari (FEM)", value=False)
    eps50 = 0.01
    if usa_py:
        eps50 = st.number_input("Deformazione eps50 per argille", 0.001, 0.05, 0.01, 0.001, format="%.4f")

    st.header("Stratigrafia")
    st.caption("Righe: spessore,gamma_dry,gamma_sat,phi,cu,k_kN_m3")
    stratigrafia_csv = st.text_area("Stratigrafia", value=str(defaults["stratigrafia_csv"]), height=150)

dati = DatiPalo(
    diametro=diametro, lunghezza=lunghezza, E=E, I=I,
    N=N, H=H,
    Nq=Nq, Nc=Nc, beta=beta, alpha=alpha, falda=falda,
    kh=kh, kv=kv,
    n_elementi=n_elementi,
    stratigrafia_csv=stratigrafia_csv
)

errors = valida_dati(dati)
if errors:
    for error in errors:
        st.error(error)
    st.stop()

try:
    results = calcola_palo(dati, usa_py=usa_py, eps50=eps50)
except Exception as e:
    st.error(f"Errore durante il calcolo: {e}")
    st.exception(e)
    st.stop()

current_input = dict(dati.__dict__)
del current_input['stratigrafia_csv'] # Keep default in UI
current_input['stratigrafia_csv'] = stratigrafia_csv

with st.sidebar:
    st.divider()
    st.download_button(
        label="Esporta input JSON",
        data=export_json(current_input),
        file_name="input_palo.json",
        mime="application/json",
        use_container_width=True,
    )

# --- Main Page ---
summary = riepilogo_progetto(dati, results)

m1, m2, m3 = st.columns(3)
m1.metric("Rapporto D/C Governanate", f"{summary['governing_ratio']:.2f}", delta=summary['governing_combination'])
m2.metric("Qamm Statica [kN]", f"{summary['Qamm_static_kN']:.0f}")
m3.metric("Spost. Max [mm]", f"{summary['y_max_mm']:.1f}")

tabs = st.tabs([
    "Sintesi e Verifiche",
    "Modello e Tensioni",
    "Risposta Laterale",
    "Confronto FD/FEM",
    "Note e Avvisi",
])

with tabs[0]:
    st.markdown('<div class="section-card">', unsafe_allow_html=True)
    st.subheader("Esiti dei Controlli Principali")
    esiti_df = valuta_esiti(dati, results)
    st.dataframe(esiti_df, use_container_width=True, hide_index=True)
    st.markdown('</div>', unsafe_allow_html=True)

    st.markdown('<div class="section-card">', unsafe_allow_html=True)
    st.subheader("Tabella Domanda/Capacità")
    st.dataframe(tabella_domanda_capacita(dati, results).round(2), use_container_width=True, hide_index=True)
    st.markdown('</div>', unsafe_allow_html=True)

with tabs[1]:
    st.markdown('<div class="section-card">', unsafe_allow_html=True)
    st.subheader("Geometria, Stratigrafia e Falda")
    st.plotly_chart(figura_geometria(dati, results), use_container_width=True)
    st.markdown('</div>', unsafe_allow_html=True)

    st.markdown('<div class="section-card">', unsafe_allow_html=True)
    st.subheader("Profilo Tensioni Efficaci")
    st.dataframe(tabella_tensioni_effettive(dati).round(2), use_container_width=True, hide_index=True)
    st.markdown('</div>', unsafe_allow_html=True)

with tabs[2]:
    st.markdown('<div class="section-card">', unsafe_allow_html=True)
    st.subheader("Diagrammi Spostamenti e Momenti (Modello FD)")
    st.plotly_chart(figura_spostamenti(results), use_container_width=True)
    st.markdown('</div>', unsafe_allow_html=True)

with tabs[3]:
    st.markdown('<div class="section-card">', unsafe_allow_html=True)
    st.subheader("Confronto Risultati Laterali: Differenze Finite vs FEM")
    st.dataframe(tabella_confronto_palo(results), use_container_width=True, hide_index=True)
    st.plotly_chart(figura_confronto_palo(results), use_container_width=True)
    if 'errore' in results['statico']['lat_fem']:
        st.error(f"Errore nel solutore FEM: {results['statico']['lat_fem']['errore']}")
    st.markdown('</div>', unsafe_allow_html=True)

with tabs[4]:
    st.markdown('<div class="section-card">', unsafe_allow_html=True)
    st.subheader("Avvisi Automatici")
    warnings = genera_warning(dati, results)
    if warnings:
        for warning in warnings:
            st.warning(warning)
    else:
        st.success("Nessun avviso automatico rilevato.")

    st.subheader("Note Tecniche")
    notes = genera_note(dati, results)
    for note in notes:
        st.info(note)
    st.markdown('</div>', unsafe_allow_html=True)