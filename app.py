# -*- coding: utf-8 -*-
import json
import streamlit as st
from src import (
    DatiPlintoPali, DEFAULT_STRAT, valida_dati, calcola_plinto_pali,
    tabella_sintesi, tabella_pali_comparativa, genera_warning, export_json,
    figura_3d_plinto_pali, figura_stratigrafia, figura_mesh_fem, genera_verifiche_df,
    figura_geometria, figura_comparativa_rigido_fem
)

DEFAULTS = {
    'B': 5.0, 'L': 5.0, 'spessore_plinto': 1.0, 'E_cls_MPa': 30000.0,
    'n_x': 3, 'n_y': 3, 'interasse_x': 1.8, 'interasse_y': 1.8,
    'diametro_palo': 0.6, 'lunghezza_palo': 20.0,
    'Nq': 25.0, 'Nc': 9.0, 'beta': 0.30, 'alpha': 0.70, 'gamma_sicurezza': 2.5,
    'N': 6000.0, 'Mx': 1500.0, 'My': 1000.0, 'kh': 0.15, 'kv': 0.05, 'falda': 2.0,
    'stratigrafia_csv': DEFAULT_STRAT,
}

st.set_page_config(page_title='PlintoPali Ultimate', layout='wide')
st.title('PlintoPali Ultimate - Analisi Strutturale & Geotecnica 3D')
st.markdown("Strumento avanzato per l'analisi di fondazioni profonde. Include stima automatica delle rigidezze, effetto gruppo, interazione terreno-struttura con graticcio **OpenSeesPy** e metodo S&T.")

with st.sidebar:
    st.header('📂 Import / Export')
    up = st.file_uploader('Carica configurazione (JSON)', type=['json'])
    defaults = DEFAULTS.copy()
    if up is not None:
        try:
            defaults.update(json.load(up))
            st.success('Dati caricati!')
        except Exception:
            st.error('File non valido.')

    st.header('📐 Geometria Plinto')
    B = st.number_input('Dimensione B [m]', 1.0, 50.0, float(defaults['B']), 0.5)
    L = st.number_input('Dimensione L [m]', 1.0, 50.0, float(defaults['L']), 0.5)
    spessore_plinto = st.number_input('Spessore Plinto H [m]', 0.3, 5.0, float(defaults['spessore_plinto']), 0.1)
    E_cls_MPa = st.number_input('Modulo E Cls [MPa]', 10000.0, 60000.0, float(defaults['E_cls_MPa']), 1000.0)

    st.header('📍 Geometria Pali')
    n_x = st.number_input('Pali X', 1, 15, int(defaults['n_x']), 1)
    n_y = st.number_input('Pali Y', 1, 15, int(defaults['n_y']), 1)
    interasse_x = st.number_input('Interasse X [m]', 0.5, 10.0, float(defaults['interasse_x']), 0.1)
    interasse_y = st.number_input('Interasse Y [m]', 0.5, 10.0, float(defaults['interasse_y']), 0.1)
    diametro_palo = st.number_input('Diametro D [m]', 0.2, 3.0, float(defaults['diametro_palo']), 0.1)
    lunghezza_palo = st.number_input('Lunghezza L [m]', 2.0, 100.0, float(defaults['lunghezza_palo']), 1.0)

    st.header('⚙️ Azioni Baricentro')
    N = st.number_input('N [kN]', 0.0, 1e7, float(defaults['N']), 500.0)
    Mx = st.number_input('Mx [kNm]', -1e6, 1e6, float(defaults['Mx']), 100.0)
    My = st.number_input('My [kNm]', -1e6, 1e6, float(defaults['My']), 100.0)

    st.header('🌍 Geotecnica')
    Nq = st.number_input('Nq (punta sabbia)', 1.0, 100.0, float(defaults['Nq']))
    Nc = st.number_input('Nc (punta argilla)', 1.0, 20.0, float(defaults['Nc']))
    beta = st.number_input('β (attrito)', 0.0, 2.0, float(defaults['beta']))
    alpha = st.number_input('α (aderenza)', 0.0, 2.0, float(defaults['alpha']))
    gamma_sicurezza = st.number_input('Fattore di Sicurezza FS', 1.0, 5.0, float(defaults['gamma_sicurezza']))
    falda = st.number_input('Profondità Falda [m]', 0.0, 100.0, float(defaults['falda']))
    
    st.subheader('Stratigrafia ed Edometrico')
    st.caption('Modifica la stratigrafia direttamente nella tabella (copia/incolla da Excel supportato).')
    
    # Parse default stratigraphy for the editor
    initial_strat_df_for_editor, _ = parse_stratigrafia(defaults['stratigrafia_csv'])
    editable_strat_df_cols = ['spessore_m', 'gamma_dry', 'gamma_sat', 'phi_deg', 'cu_kPa', 'E_ed_kPa']
    
    edited_strat_df = st.data_editor(
        initial_strat_df_for_editor[editable_strat_df_cols],
        num_rows="dynamic",
        use_container_width=True,
        column_config={
            "spessore_m": st.column_config.NumberColumn("Spessore [m]", format="%.1f"),
            "gamma_dry": st.column_config.NumberColumn("γ_dry [kN/m³]", format="%.1f"),
            "gamma_sat": st.column_config.NumberColumn("γ_sat [kN/m³]", format="%.1f"),
            "phi_deg": st.column_config.NumberColumn("φ [°]", format="%.0f"),
            "cu_kPa": st.column_config.NumberColumn("cᵤ [kPa]", format="%.0f"),
            "E_ed_kPa": st.column_config.NumberColumn("E_ed [kPa]", format="%.0f"),
        },
        key="strat_editor"
    )
    stratigrafia_csv = edited_strat_df.to_csv(header=False, index=False, float_format='%.1f')

# Costruzione istanza
d = DatiPlintoPali(
    B, L, spessore_plinto, E_cls_MPa, int(n_x), int(n_y), interasse_x, interasse_y,
    diametro_palo, lunghezza_palo, Nq, Nc, beta, alpha, gamma_sicurezza, N, Mx, My,
    0.0, 0.0, falda, stratigrafia_csv
)

err = valida_dati(d)
if err:
    for e in err: st.error(e)
    st.stop()

try:
    with st.spinner("Analisi FEM e geotecnica in corso..."):
        r = calcola_plinto_pali(d)
    
    verifiche_df = genera_verifiche_df(d, r)
    warnings_count = len(verifiche_df[verifiche_df['Esito'].isin(["NON VERIFICATO", "ATTENZIONE"])])
    current_input = {k: v for k, v in d.__dict__.items()}
    
    # Intestazione e Avvisi Rapidi
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Pali Totali", f"{d.n_x * d.n_y}")
    c2.metric("Portata Amm. (Gruppo)", f"{r['Qamm_effettiva_palo']:.0f} kN")
    c3.metric("Cedimento Stimato", f"{r['cedimento_gruppo_mm']:.2f} mm") # Ensure 2 decimal places
    c4.metric("Stato Verifiche", f"{warnings_count} ⚠" if warnings_count > 0 else "OK ✓")

            
    # Gestione Tabs
    t1, t2, t3, t4, t5, t6 = st.tabs([
        '🏗️ Vista 3D & Terreno', 
        '📊 Confronto Reazioni', 
        '🛠️ Sollecitazioni Plinto', 
        '📑 Report Tabellare',
        '💾 Salva Dati'
    ])
    
    # Importa il modulo di reporting solo se necessario
    try:
        from reporting import crea_report_word
        reporting_enabled = True # type: ignore
    except ImportError:
        reporting_enabled = False

    with t1:
        st.subheader("Geometria 3D ed Esplorazione Geotecnica")
        col_3d, col_strat = st.columns([2, 1])
        with col_3d:
            st.plotly_chart(figura_3d_plinto_pali(d), use_container_width=True)
        with col_strat:
            st.plotly_chart(figura_stratigrafia(r), use_container_width=True) # Ensure 2 decimal places in plot
            eff = r['efficienza_gruppo']
            st.info(f"**Mod. Assiale Kv:** {r['k_v_calcolato']:.0f} kN/m\n\n**Interazione:** {eff['stato']} ($\eta$ = {eff['eta']:.2f})")
            
    with t2:
        st.subheader("Analisi delle Reazioni: Plinto Infinitamente Rigido vs Flessibile")
        st.plotly_chart(figura_comparativa_rigido_fem(r), use_container_width=True)
        st.plotly_chart(figura_geometria(d, r), use_container_width=True)

    with t3:
        st.header("Verifiche Strutturali della Piastra di Fondazione")
        
        st.subheader("1. Modello Flessibile FEM (Graticcio OpenSeesPy)")
        col_mesh, col_data = st.columns([2, 1])
        with col_mesh:
            st.plotly_chart(figura_mesh_fem(d, r), use_container_width=True)
        with col_data:
            m_max = max(r['statico_fem']['M_radice_kNm'])
            st.metric("Momento Max Graticcio (Incastro)", f"{m_max:.1f} kNm")
            st.markdown("Mostra i cedimenti puntuali ai nodi in mm e calcola la deformata elastica degli elementi trave.")
            
        st.divider()
        
        st.subheader("2. Modello a Puntoni e Tiranti (Ipotesi Strut-and-Tie / Rigido)")
        st_tie = r['strut_and_tie']
        st.markdown(f"""
        Assumendo l'altezza utile $d = {st_tie['d_utile_m']:.2f}$ m:
        * Trazione Totale Inferiore X ($T_x$): **{st_tie['Tx_kN']:.0f} kN**
        * Trazione Totale Inferiore Y ($T_y$): **{st_tie['Ty_kN']:.0f} kN**
        """)

    with t4:
        st.subheader("Tabelle Dettagliate per Post-Processing")
        df_comp = tabella_pali_comparativa(r)
        st.dataframe(df_comp, use_container_width=True)
        st.download_button(
            'Scarica Tabella (CSV)', 
            df_comp.to_csv(index=False, float_format='%.2f').encode('utf-8'), # Ensure 2 decimal places in CSV
            'plintopali_reazioni.csv', 
            'text/csv'
        )
    
    with t4: # New Verifiche tab
        st.subheader("Verifiche di Progetto")
        st.dataframe(verifiche_df, use_container_width=True)

    with t5:
        st.subheader("Salvataggio Impostazioni")
        st.download_button(
            'Scarica Input in Formato JSON', 
            export_json(current_input), 
            'plintopali_input.json', 
            'application/json'
        )
        
        if reporting_enabled:
            st.divider()
            st.subheader("Generazione Relazione Tecnica")
            if st.button("Genera Relazione (.docx)"):
                with st.spinner("Creazione del documento Word..."):
                    report_bytes = crea_report_word(d, r, df_comp, verifiche_df)
                    st.download_button(
                        label="Scarica Relazione Word",
                        data=report_bytes,
                        file_name="Relazione_PlintoPali.docx",
                        mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document"
                    )

except Exception as e:
    st.error(f"Errore critico durante la risoluzione del modello: {e}")