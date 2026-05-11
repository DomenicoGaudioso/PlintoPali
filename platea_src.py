# -*- coding: utf-8 -*-
from __future__ import annotations
from dataclasses import dataclass
from typing import List, Dict
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import openseespy.opensees as ops

@dataclass(frozen=True)
class DatiPlatea:
    B: float
    L: float
    spessore: float
    E_cls_MPa: float
    k_winkler_kPa_m: float
    mesh_size: float
    pilastri_df: pd.DataFrame
    q_distribuito_kPa: float = 0.0
    poisson: float = 0.2

def valida_dati_platea(d: DatiPlatea) -> List[str]:
    err = []
    if d.B <= 0 or d.L <= 0 or d.spessore <= 0:
        err.append('Le dimensioni della platea devono essere positive.')
    if d.E_cls_MPa <= 0:
        err.append('Il modulo elastico del calcestruzzo deve essere positivo.')
    if d.k_winkler_kPa_m <= 0:
        err.append('Il modulo di Winkler deve essere positivo.')
    if d.mesh_size <= 0:
        err.append('La dimensione della mesh deve essere positiva.')
    if 'x' not in d.pilastri_df.columns or 'y' not in d.pilastri_df.columns:
        err.append("La tabella pilastri deve contenere le colonne 'x' e 'y'.")
    else:
        if not d.pilastri_df.empty and (d.pilastri_df['x'].max() > d.B or d.pilastri_df['y'].max() > d.L):
            err.append('Posizione pilastro fuori dalla geometria della platea.')
    return err

def calcola_platea_fem(d: DatiPlatea) -> Dict:
    """
    Analisi FEM di una platea su suolo elastico (Winkler) con OpenSeesPy.
    La platea è modellata con elementi ShellMITC4.
    Il terreno è modellato con molle verticali indipendenti.
    """
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)

    # --- 1. Creazione Mesh e Nodi ---
    nx = int(np.ceil(d.B / d.mesh_size)) + 1
    ny = int(np.ceil(d.L / d.mesh_size)) + 1
    x_coords = np.linspace(0, d.B, nx)
    y_coords = np.linspace(0, d.L, ny)
    
    node_tags = np.arange(1, nx * ny + 1).reshape((ny, nx))
    
    for i in range(ny):
        for j in range(nx):
            node_tag = int(node_tags[i, j])
            ops.node(node_tag, x_coords[j], y_coords[i], 0.0)

    # --- 2. Materiali ---
    E_cls_kPa = d.E_cls_MPa * 1000
    nu = d.poisson
    # Materiale elastico isotropo per la platea
    mat_tag = 1
    ops.nDMaterial('ElasticIsotropic', mat_tag, E_cls_kPa, nu)

    # --- 3. Elementi Shell ---
    ele_type = 'ShellMITC4'
    sect_tag = 1
    ops.section('PlateFiber', sect_tag, mat_tag, d.spessore)

    element_tags = np.arange(1, (nx - 1) * (ny - 1) + 1).reshape((ny - 1, nx - 1))
    for i in range(ny - 1):
        for j in range(nx - 1):
            ele_tag = int(element_tags[i, j])
            n1 = int(node_tags[i, j])
            n2 = int(node_tags[i, j+1])
            n3 = int(node_tags[i+1, j+1])
            n4 = int(node_tags[i+1, j])
            ops.element(ele_type, ele_tag, n1, n2, n3, n4, sect_tag)

    # --- 4. Vincoli (Molle di Winkler) ---
    # Calcolo area di influenza per ogni nodo
    node_area = d.mesh_size * d.mesh_size 
    k_node = d.k_winkler_kPa_m * node_area

    # Materiale elastico per le molle
    spring_mat_tag_start = 1000
    for i in range(ny):
        for j in range(nx):
            node_tag = int(node_tags[i, j])
            # Aggiusta l'area per i nodi sui bordi e angoli
            area_multiplier = 1.0
            if (i==0 or i==ny-1) and (j==0 or j==nx-1): # Angolo
                area_multiplier = 0.25
            elif i==0 or i==ny-1 or j==0 or j==nx-1: # Bordo
                area_multiplier = 0.5
            
            k_node_eff = d.k_winkler_kPa_m * (d.mesh_size**2) * area_multiplier

            spring_mat_tag = spring_mat_tag_start + node_tag
            ops.uniaxialMaterial('Elastic', spring_mat_tag, k_node_eff)
            
            # Elemento molla (zeroLength)
            spring_ele_tag = spring_mat_tag # Stesso tag per semplicità
            ops.element('zeroLength', spring_ele_tag, node_tag, '-mat', spring_mat_tag, '-dir', 3)
            # Fissa l'altro capo della molla
            ops.fix(node_tag, 0, 0, 1, 0, 0, 0) # Vincola solo la traslazione Z

    # --- 5. Carichi ---
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)

    for _, pilastro in d.pilastri_df.iterrows():
        px, py = pilastro['x'], pilastro['y']
        # Trova il nodo della mesh più vicino
        idx_x = (np.abs(x_coords - px)).argmin()
        idx_y = (np.abs(y_coords - py)).argmin()
        target_node = int(node_tags[idx_y, idx_x])
        
        P = -float(pilastro.get('P_kN', 0.0)) # Negativo perché Z è verso l'alto
        Mx = float(pilastro.get('Mx_kNm', 0.0))
        My = float(pilastro.get('My_kNm', 0.0))
        ops.load(target_node, 0.0, 0.0, P, Mx, My, 0.0)

    # Carico distribuito
    if d.q_distribuito_kPa != 0.0:
        for i in range(ny):
            for j in range(nx):
                node_tag = int(node_tags[i, j])
                
                # Calcolo area di influenza per il carico distribuito
                area_multiplier = 1.0
                if (i==0 or i==ny-1) and (j==0 or j==nx-1): # Angolo
                    area_multiplier = 0.25
                elif i==0 or i==ny-1 or j==0 or j==nx-1: # Bordo
                    area_multiplier = 0.5
                
                node_area = (d.mesh_size**2) * area_multiplier
                nodal_force = -d.q_distribuito_kPa * node_area # Negativo perché Z è verso l'alto
                ops.load(node_tag, 0.0, 0.0, nodal_force, 0.0, 0.0, 0.0)

    # --- 6. Analisi ---
    ops.system('BandGeneral')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')

    if ops.analyze(1) != 0:
        raise RuntimeError("Analisi OpenSees fallita. Controllare carichi e vincoli.")

    # --- 7. Estrazione Risultati ---
    cedimenti = np.zeros((ny, nx))
    pressioni = np.zeros((ny, nx))
    Mxx = np.zeros((ny - 1, nx - 1))
    Myy = np.zeros((ny - 1, nx - 1))
    Mxy = np.zeros((ny - 1, nx - 1))

    for i in range(ny):
        for j in range(nx):
            node_tag = int(node_tags[i, j])
            cedimenti[i, j] = ops.nodeDisp(node_tag, 3) * 1000 # in mm
            
            # Calcolo pressione dal cedimento
            area_multiplier = 1.0
            if (i==0 or i==ny-1) and (j==0 or j==nx-1): area_multiplier = 0.25
            elif i==0 or i==ny-1 or j==0 or j==nx-1: area_multiplier = 0.5
            k_node_eff = d.k_winkler_kPa_m * (d.mesh_size**2) * area_multiplier
            force = k_node_eff * abs(ops.nodeDisp(node_tag, 3))
            pressioni[i,j] = force / (d.mesh_size**2 * area_multiplier)

    for i in range(ny - 1):
        for j in range(nx - 1):
            ele_tag = int(element_tags[i, j])
            forces = ops.eleResponse(ele_tag, 'force')
            # I momenti sono per unità di lunghezza
            Mxx[i, j] = forces[0]
            Myy[i, j] = forces[1]
            Mxy[i, j] = forces[2]

    return {
        'x_coords': x_coords,
        'y_coords': y_coords,
        'node_tags': node_tags,
        'element_tags': element_tags,
        'cedimenti_mm': cedimenti,
        'pressioni_kPa': pressioni,
        'Mxx_kNm_m': Mxx,
        'Myy_kNm_m': Myy,
        'Mxy_kNm_m': Mxy,
    }

def figura_geometria_platea(d: DatiPlatea, r: Dict) -> go.Figure:
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[0, d.B, d.B, 0, 0], y=[0, 0, d.L, d.L, 0], fill='toself', mode='lines', name='Platea'))
    fig.add_trace(go.Scatter(x=d.pilastri_df['x'], y=d.pilastri_df['y'], mode='markers+text', name='Pilastri',
                             marker=dict(size=10, color='red'), text=d.pilastri_df['P_kN'].astype(str) + ' kN', textposition='top center'))
    
    # Aggiungi mesh
    for i in range(len(r['x_coords'])):
        fig.add_shape(type='line', x0=r['x_coords'][i], y0=0, x1=r['x_coords'][i], y1=d.L, line=dict(color='lightgrey', width=1))
    for i in range(len(r['y_coords'])):
        fig.add_shape(type='line', x0=0, y0=r['y_coords'][i], x1=d.B, y1=r['y_coords'][i], line=dict(color='lightgrey', width=1))

    fig.update_layout(title='Geometria e Mesh', yaxis_scaleanchor="x", xaxis_constrain='domain')
    return fig

def figura_risultati_platea(d: DatiPlatea, r: Dict, z_key: str, title: str) -> go.Figure:
    is_moment = 'M' in z_key
    x = r['x_coords']
    y = r['y_coords']
    z = r[z_key]

    if is_moment:
        # I momenti sono al centro dell'elemento, quindi le coordinate sono shiftate
        x = (r['x_coords'][:-1] + r['x_coords'][1:]) / 2
        y = (r['y_coords'][:-1] + r['y_coords'][1:]) / 2

    fig = go.Figure(data=go.Contour(
        z=z,
        x=x,
        y=y,
        colorscale='Viridis',
        colorbar_title=title,
        contours_coloring='lines',
        line_width=2,
    ))
    fig.add_trace(go.Contour(z=z, x=x, y=y, colorscale='Viridis', showscale=False, contours_coloring='heatmap'))
    fig.update_layout(title=title, yaxis_scaleanchor="x", xaxis_constrain='domain')
    return fig