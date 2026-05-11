# -*- coding: utf-8 -*-
import io
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from docx import Document
from docx.shared import Inches, Pt


def add_df_to_doc(doc, df: pd.DataFrame, font_size: int = 9):
    """Aggiunge un DataFrame pandas a un documento docx."""
    df = df.round(2)
    t = doc.add_table(df.shape[0] + 1, df.shape[1])
    t.style = 'Table Grid'
    for j, col_name in enumerate(df.columns):
        t.cell(0, j).text = col_name
    for i, row in enumerate(df.itertuples(), start=1):
        for j, val in enumerate(row[1:]):
            cell = t.cell(i, j)
            cell.text = str(val)
            for paragraph in cell.paragraphs:
                for run in paragraph.runs:
                    run.font.size = Pt(font_size)


def plot_risultati_platea_mpl(risultati, z_key, title):
    """Genera un grafico a contorni con Matplotlib."""
    fig, ax = plt.subplots(figsize=(8, 7))

    is_moment = 'M' in z_key
    x = risultati['x_coords']
    y = risultati['y_coords']
    z = risultati[z_key]

    if is_moment:
        x = (risultati['x_coords'][:-1] + risultati['x_coords'][1:]) / 2
        y = (risultati['y_coords'][:-1] + risultati['y_coords'][1:]) / 2

    # Matplotlib's contourf expects Z with shape (Y.shape[0], X.shape[0])
    # I dati da OpenSees sono (ny, nx), che corrisponde a (len(y), len(x))
    # quindi non serve la trasposizione se le coordinate sono corrette.
    contour = ax.contourf(x, y, z, cmap='viridis', levels=15)
    cbar = fig.colorbar(contour, ax=ax)
    
    try:
        cbar.set_label(title.split('[')[1].split(']')[0])
    except IndexError:
        cbar.set_label(title)

    ax.contour(x, y, z, colors='k', linewidths=0.5, levels=15)

    ax.set_title(title)
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_aspect('equal', adjustable='box')
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150)
    plt.close(fig)
    buf.seek(0)
    return buf


def crea_report_word_platea(dati, risultati):
    """Crea un documento Word con i risultati dell'analisi della platea."""
    doc = Document()
    doc.add_heading('Relazione di Calcolo - Platea di Fondazione (FEM)', 0)
    doc.add_paragraph(f"Software: PlateaFEM")
    doc.add_paragraph(f"Data: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    doc.add_heading('1. Dati di Input', level=1)
    doc.add_paragraph(f"Dimensioni platea: {dati.B:.2f}m (X) x {dati.L:.2f}m (Y)")
    doc.add_paragraph(f"Spessore: {dati.spessore:.2f} m")
    doc.add_paragraph(f"Modulo E calcestruzzo: {dati.E_cls_MPa:.0f} MPa")
    doc.add_paragraph(f"Modulo di Winkler: {dati.k_winkler_kPa_m:.0f} kPa/m")
    doc.add_paragraph(f"Dimensione mesh: {dati.mesh_size:.2f} m")
    if dati.q_distribuito_kPa != 0:
        doc.add_paragraph(f"Carico distribuito: {dati.q_distribuito_kPa:.2f} kPa")

    if not dati.pilastri_df.empty:
        doc.add_heading('Carichi sui Pilastri', level=2)
        add_df_to_doc(doc, dati.pilastri_df)

    doc.add_heading('2. Risultati di Sintesi', level=1)
    cedimento_max = risultati['cedimenti_mm'].max()
    cedimento_min = risultati['cedimenti_mm'].min()
    pressione_max = risultati['pressioni_kPa'].max()
    mxx_max = risultati['Mxx_kNm_m'].max()
    myy_max = risultati['Myy_kNm_m'].max()

    doc.add_paragraph(f"Cedimento massimo: {cedimento_max:.2f} mm")
    doc.add_paragraph(f"Cedimento minimo: {cedimento_min:.2f} mm")
    doc.add_paragraph(f"Pressione massima sul terreno: {pressione_max:.1f} kPa")
    doc.add_paragraph(f"Momento Mxx massimo: {mxx_max:.1f} kNm/m")
    doc.add_paragraph(f"Momento Myy massimo: {myy_max:.1f} kNm/m")

    doc.add_heading('3. Mappe dei Risultati', level=1)

    doc.add_paragraph("Mappa dei cedimenti:")
    img_ced = plot_risultati_platea_mpl(risultati, 'cedimenti_mm', 'Cedimenti [mm]')
    doc.add_picture(img_ced, width=Inches(6.0))

    doc.add_paragraph("Mappa delle pressioni sul terreno:")
    img_press = plot_risultati_platea_mpl(risultati, 'pressioni_kPa', 'Pressioni [kPa]')
    doc.add_picture(img_press, width=Inches(6.0))

    doc.add_paragraph("Mappa dei momenti flettenti Mxx:")
    img_mxx = plot_risultati_platea_mpl(risultati, 'Mxx_kNm_m', 'Momento Mxx [kNm/m]')
    doc.add_picture(img_mxx, width=Inches(6.0))

    doc.add_paragraph("Mappa dei momenti flettenti Myy:")
    img_myy = plot_risultati_platea_mpl(risultati, 'Myy_kNm_m', 'Momento Myy [kNm/m]')
    doc.add_picture(img_myy, width=Inches(6.0))

    doc_io = io.BytesIO()
    doc.save(doc_io)
    doc_io.seek(0)
    return doc_io.getvalue()