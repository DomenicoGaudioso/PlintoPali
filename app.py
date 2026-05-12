# -*- coding: utf-8 -*-
from __future__ import annotations

from dataclasses import asdict
import json

import pandas as pd
import streamlit as st

from report import create_word_report
from src import (
    DEFAULT_STRAT,
    DatiPlintoPali,
    calcola_plinto_pali,
    export_json,
    figura_3d_plinto_pali,
    figura_comparativa,
    figura_geometria,
    figura_mesh_fem,
    figura_output,
    figura_stratigrafia,
    genera_note,
    genera_verifiche_df,
    genera_warning,
    tabella_pali_comparativa,
    tabella_sintesi,
    valida_dati,
)


DEFAULTS = {
    "B": 5.0,
    "L": 5.0,
    "spessore_plinto": 1.0,
    "E_cls_MPa": 30000.0,
    "n_x": 2,
    "n_y": 2,
    "interasse_x": 3.0,
    "interasse_y": 3.0,
    "diametro_palo": 0.6,
    "lunghezza_palo": 12.0,
    "Nq": 35.0,
    "Nc": 9.0,
    "beta": 0.35,
    "alpha": 0.7,
    "gamma_sicurezza": 2.5,
    "N": 2400.0,
    "Mx": 250.0,
    "My": 180.0,
    "kh": 0.06,
    "kv": 0.0,
    "falda": 100.0,
    "stratigrafia_csv": DEFAULT_STRAT,
    "layout_pali": None,
}


def _grid_layout_rows(n_x: int, n_y: int, interasse_x: float, interasse_y: float) -> list[dict[str, float]]:
    rows = []
    x0 = -(n_x - 1) * interasse_x / 2.0
    y0 = -(n_y - 1) * interasse_y / 2.0
    for iy in range(n_y):
        for ix in range(n_x):
            rows.append({"x": x0 + ix * interasse_x, "y": y0 + iy * interasse_y})
    return rows


def _load_defaults(uploaded) -> dict:
    values = DEFAULTS.copy()
    if uploaded is None:
        return values
    payload = json.loads(uploaded.getvalue().decode("utf-8"))
    raw = payload.get("input", payload)
    for key in DEFAULTS:
        if key in raw:
            values[key] = raw[key]
    return values


st.set_page_config(page_title="PlintoPali", layout="wide")
st.title("PlintoPali")
st.caption("Analisi statica e sismica di plinti su pali con confronto rigido/FEM.")

with st.sidebar:
    st.header("Import / Export input")
    uploaded = st.file_uploader("Reimporta input JSON", type=["json"])
    try:
        defaults = _load_defaults(uploaded)
        if uploaded is not None:
            st.success("Input importati.")
    except Exception as exc:
        defaults = DEFAULTS.copy()
        st.error(f"JSON non valido: {exc}")

    st.header("Geometria plinto")
    B = st.number_input("B [m]", 0.5, 50.0, float(defaults["B"]), 0.1)
    L = st.number_input("L [m]", 0.5, 50.0, float(defaults["L"]), 0.1)
    spessore_plinto = st.number_input("Spessore plinto [m]", 0.2, 5.0, float(defaults["spessore_plinto"]), 0.05)
    E_cls_MPa = st.number_input("E cls [MPa]", 10000.0, 60000.0, float(defaults["E_cls_MPa"]), 500.0)

    st.header("Pali")
    n_x = st.number_input("Numero pali in X [-]", 1, 10, int(defaults["n_x"]), 1)
    n_y = st.number_input("Numero pali in Y [-]", 1, 10, int(defaults["n_y"]), 1)
    interasse_x = st.number_input("Interasse X [m]", 0.5, 20.0, float(defaults["interasse_x"]), 0.1)
    interasse_y = st.number_input("Interasse Y [m]", 0.5, 20.0, float(defaults["interasse_y"]), 0.1)
    diametro_palo = st.number_input("Diametro palo [m]", 0.2, 3.0, float(defaults["diametro_palo"]), 0.05)
    lunghezza_palo = st.number_input("Lunghezza palo [m]", 1.0, 80.0, float(defaults["lunghezza_palo"]), 0.5)
    layout_importato = defaults.get("layout_pali") or []
    usa_layout_custom = st.checkbox("Usa layout tabellare/irregolare", value=bool(layout_importato))
    if usa_layout_custom:
        layout_base = layout_importato or _grid_layout_rows(int(n_x), int(n_y), float(interasse_x), float(interasse_y))
        layout_df = pd.DataFrame(layout_base)
        if layout_df.empty:
            layout_df = pd.DataFrame(columns=["x", "y"])
        layout_df = layout_df.reindex(columns=["x", "y"])
        layout_editato = st.data_editor(
            layout_df,
            num_rows="dynamic",
            use_container_width=True,
            hide_index=True,
            column_config={
                "x": st.column_config.NumberColumn("x [m]", step=0.05, format="%.3f"),
                "y": st.column_config.NumberColumn("y [m]", step=0.05, format="%.3f"),
            },
        )
        layout_pali = [
            {"x": float(row["x"]), "y": float(row["y"])}
            for _, row in layout_editato.dropna(subset=["x", "y"]).iterrows()
        ]
    else:
        layout_pali = None

    st.header("Capacita geotecnica")
    Nq = st.number_input("Nq [-]", 1.0, 100.0, float(defaults["Nq"]), 1.0)
    Nc = st.number_input("Nc [-]", 1.0, 20.0, float(defaults["Nc"]), 0.5)
    beta = st.number_input("beta [-]", 0.0, 2.0, float(defaults["beta"]), 0.05)
    alpha = st.number_input("alpha [-]", 0.0, 2.0, float(defaults["alpha"]), 0.05)
    gamma_sicurezza = st.number_input("Fattore sicurezza [-]", 1.01, 5.0, float(defaults["gamma_sicurezza"]), 0.05)

    st.header("Azioni")
    N = st.number_input("N [kN]", 0.0, 100000.0, float(defaults["N"]), 50.0)
    Mx = st.number_input("Mx [kNm]", -50000.0, 50000.0, float(defaults["Mx"]), 10.0)
    My = st.number_input("My [kNm]", -50000.0, 50000.0, float(defaults["My"]), 10.0)
    kh = st.number_input("kh [-]", 0.0, 1.0, float(defaults["kh"]), 0.01)
    kv = st.number_input("kv [-]", 0.0, 1.0, float(defaults["kv"]), 0.01)
    falda = st.number_input("Profondita falda [m]", 0.0, 100.0, float(defaults["falda"]), 0.1)

    st.header("Stratigrafia")
    stratigrafia_csv = st.text_area(
        "Righe: spessore,gamma_dry,gamma_sat,phi,cu,E_ed",
        value=str(defaults["stratigrafia_csv"]),
        height=130,
    )

dati = DatiPlintoPali(
    B=B,
    L=L,
    spessore_plinto=spessore_plinto,
    E_cls_MPa=E_cls_MPa,
    n_x=int(n_x),
    n_y=int(n_y),
    interasse_x=interasse_x,
    interasse_y=interasse_y,
    diametro_palo=diametro_palo,
    lunghezza_palo=lunghezza_palo,
    Nq=Nq,
    Nc=Nc,
    beta=beta,
    alpha=alpha,
    gamma_sicurezza=gamma_sicurezza,
    N=N,
    Mx=Mx,
    My=My,
    kh=kh,
    kv=kv,
    falda=falda,
    stratigrafia_csv=stratigrafia_csv,
    layout_pali=layout_pali,
)

errors = valida_dati(dati)
if errors:
    for error in errors:
        st.error(error)
    st.stop()

risultati = calcola_plinto_pali(dati)
sintesi = tabella_sintesi(risultati)
confronto = tabella_pali_comparativa(risultati)
verifiche = genera_verifiche_df(dati, risultati)

with st.sidebar:
    st.divider()
    st.download_button(
        "Salva input JSON",
        data=export_json(asdict(dati)),
        file_name="plinto_pali_input.json",
        mime="application/json",
    )
    try:
        st.download_button(
            "Scarica relazione Word",
            data=create_word_report(dati, risultati),
            file_name="relazione_plinto_pali.docx",
            mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
        )
    except Exception as exc:
        st.warning(f"Relazione Word non disponibile: {exc}")

cols = st.columns(4)
cols[0].metric("Pali totali", f"{len(risultati['statico']['x'])}")
cols[1].metric("Q amm palo", f"{risultati['Qamm_effettiva_palo']:.2f} kN")
cols[2].metric("FS minimo rigido", f"{float(risultati['statico']['FS'].min()):.2f}")
cols[3].metric("Cedimento gruppo", f"{risultati['cedimento_gruppo_mm']:.2f} mm")

tabs = st.tabs(["Sintesi", "Geometria", "Reazioni", "Verifiche", "Log tecnico"])

with tabs[0]:
    st.subheader("Sintesi")
    st.dataframe(sintesi, use_container_width=True, hide_index=True)
    st.subheader("Confronto pali")
    st.dataframe(confronto, use_container_width=True, hide_index=True)

with tabs[1]:
    left, right = st.columns([1.3, 1.0])
    with left:
        st.plotly_chart(figura_geometria(dati, risultati), use_container_width=True)
    with right:
        st.plotly_chart(figura_3d_plinto_pali(dati), use_container_width=True)
    st.plotly_chart(figura_mesh_fem(dati, risultati), use_container_width=True)
    st.plotly_chart(figura_stratigrafia(risultati), use_container_width=True)

with tabs[2]:
    st.plotly_chart(figura_output(risultati, "statico"), use_container_width=True)
    st.plotly_chart(figura_output(risultati, "sismico"), use_container_width=True)
    st.plotly_chart(figura_comparativa(risultati), use_container_width=True)

with tabs[3]:
    st.subheader("Verifiche")
    st.dataframe(verifiche, use_container_width=True, hide_index=True)
    if (verifiche["Esito"] == "NON VERIFICATO").any():
        st.error("Sono presenti verifiche non soddisfatte.")
    else:
        st.success("Le verifiche principali risultano soddisfatte.")

with tabs[4]:
    warnings = genera_warning(dati, risultati)
    if warnings:
        for warning in warnings:
            st.warning(warning)
    else:
        st.success("Nessuna criticita automatica rilevata.")
    for note in genera_note(dati, risultati):
        st.info(note)
