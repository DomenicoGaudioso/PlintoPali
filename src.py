# -*- coding: utf-8 -*-
from __future__ import annotations
from dataclasses import asdict
from typing import Tuple, List
import json
import math
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dataclasses import dataclass
from typing import Dict, List

GAMMA_W = 9.81
DEFAULT_STRAT = """2.0,18,20,30,0,25000
3.0,19,21,34,0,40000
5.0,20,20,0,120,60000
"""


def parse_stratigrafia(csv_text: str) -> Tuple[pd.DataFrame, List[str]]:
    """Righe: spessore,gamma_dry,gamma_sat,phi_deg,cu_kPa,k_kN_m3"""
    err, rows = [], []
    lines = [ln.strip() for ln in csv_text.splitlines() if ln.strip()]
    if not lines:
        return pd.DataFrame(columns=['spessore_m', 'gamma_dry', 'gamma_sat', 'phi_deg', 'cu_kPa', 'k_kN_m3']), ['Inserire almeno uno strato.']
    for i, line in enumerate(lines, start=1):
        parts = [p.strip() for p in line.replace(';', ',').split(',') if p.strip()]
        if len(parts) != 6:
            err.append(f'Riga {i}: usare 6 campi = spessore,gamma_dry,gamma_sat,phi,cu,k.')
            continue
        try:
            h, gd, gs, phi, cu, k = map(float, parts)
            rows.append({'spessore_m': h, 'gamma_dry': gd, 'gamma_sat': gs, 'phi_deg': phi, 'cu_kPa': cu, 'k_kN_m3': k})
        except ValueError:
            err.append(f'Riga {i}: valori non numerici.')
    df = pd.DataFrame(rows)
    if df.empty:
        return df, err or ['Stratigrafia non valida.']
    if (df['spessore_m'] <= 0).any():
        err.append('Tutti gli spessori devono essere positivi.')
    df['z_top_m'] = df['spessore_m'].cumsum() - df['spessore_m']
    df['z_bot_m'] = df['spessore_m'].cumsum()
    return df, err


def layer_at_depth(df: pd.DataFrame, z: float) -> pd.Series:
    sel = df[(df['z_top_m'] <= z) & (df['z_bot_m'] >= z)]
    if sel.empty:
        return df.iloc[-1]
    return sel.iloc[0]


def gamma_eff(layer: pd.Series, z_mid: float, falda_depth: float) -> float:
    if z_mid <= falda_depth:
        return float(layer['gamma_dry'])
    return max(float(layer['gamma_sat']) - GAMMA_W, 1.0)


def sigma_v_eff(df: pd.DataFrame, z: float, falda_depth: float) -> float:
    s = 0.0
    for _, r in df.iterrows():
        a = max(0.0, float(r['z_top_m']))
        b = min(z, float(r['z_bot_m']))
        if b <= a:
            continue
        if falda_depth <= a:
            s += (b - a) * max(float(r['gamma_sat']) - GAMMA_W, 1.0)
        elif falda_depth >= b:
            s += (b - a) * float(r['gamma_dry'])
        else:
            s += (falda_depth - a) * float(r['gamma_dry'])
            s += (b - falda_depth) * max(float(r['gamma_sat']) - GAMMA_W, 1.0)
    return s


def u_hydro(z: float, falda_depth: float) -> float:
    return GAMMA_W * max(z - falda_depth, 0.0)


def export_json(data: dict) -> bytes:
    return json.dumps(data, ensure_ascii=False, indent=2).encode('utf-8')


@dataclass(frozen=True)
class DatiPlintoPali:
    B: float
    L: float
    n_x: int
    n_y: int
    interasse_x: float
    interasse_y: float
    diametro_palo: float
    lunghezza_palo: float
    Nq: float
    Nc: float
    beta: float
    alpha: float
    gamma_sicurezza: float
    N: float
    Mx: float
    My: float
    kh: float
    kv: float
    falda: float
    stratigrafia_csv: str = DEFAULT_STRAT


def valida_dati(d: DatiPlintoPali) -> List[str]:
    err = []
    if d.B <= 0 or d.L <= 0:
        err.append('Le dimensioni del plinto devono essere positive.')
    if d.n_x <= 0 or d.n_y <= 0:
        err.append('Il numero di pali per direzione deve essere positivo.')
    if d.interasse_x <= 0 or d.interasse_y <= 0:
        err.append('Gli interassi devono essere positivi.')
    if d.lunghezza_palo <= 0:
        err.append('La lunghezza del palo deve essere positiva.')
    if d.gamma_sicurezza <= 1.0:
        err.append('La sicurezza del palo deve essere > 1.')
    _, e2 = parse_stratigrafia(d.stratigrafia_csv)
    err.extend(e2)
    return err


def coordinate_pali(d: DatiPlintoPali):
    xs = np.linspace(-(d.n_x - 1) * d.interasse_x / 2, (d.n_x - 1) * d.interasse_x / 2, d.n_x)
    ys = np.linspace(-(d.n_y - 1) * d.interasse_y / 2, (d.n_y - 1) * d.interasse_y / 2, d.n_y)
    X, Y = np.meshgrid(xs, ys)
    return X.flatten(), Y.flatten()


def capacita_singolo_palo(d: DatiPlintoPali) -> Dict[str, float]:
    df, _ = parse_stratigrafia(d.stratigrafia_csv)
    A_base = math.pi * d.diametro_palo ** 2 / 4.0
    P = math.pi * d.diametro_palo
    Qs = 0.0
    for _, r in df.iterrows():
        a = max(0.0, float(r['z_top_m']))
        b = min(d.lunghezza_palo, float(r['z_bot_m']))
        if b <= a:
            continue
        dz = b - a
        zmid = 0.5 * (a + b)
        phi = float(r['phi_deg'])
        cu = float(r['cu_kPa'])
        if phi > 0:
            qs = d.beta * sigma_v_eff(df, zmid, d.falda)
        else:
            qs = d.alpha * cu
        Qs += qs * P * dz
    tip = layer_at_depth(df, d.lunghezza_palo)
    if float(tip['phi_deg']) > 0:
        Qb = d.Nq * sigma_v_eff(df, d.lunghezza_palo, d.falda) * A_base
    else:
        Qb = d.Nc * float(tip['cu_kPa']) * A_base
    Qult = Qb + Qs
    return {'stratigrafia': df, 'Qult': Qult, 'Qamm': Qult / d.gamma_sicurezza, 'Qb': Qb, 'Qs': Qs}


def reaction_case(d: DatiPlintoPali, seismic=False):
    x, y = coordinate_pali(d)
    n = len(x)
    A = np.vstack([np.ones(n), y, x]).T
    N = d.N * (1.0 - (d.kv if seismic else 0.0))
    Mx = d.Mx * (1.0 + (d.kh if seismic else 0.0))
    My = d.My * (1.0 + (d.kh if seismic else 0.0))
    abc = np.linalg.solve(A.T @ A, np.array([N, Mx, My]))
    R = A @ abc
    return {'x': x, 'y': y, 'R': R}


def calcola_plinto_pali(d: DatiPlintoPali) -> Dict[str, object]:
    cap = capacita_singolo_palo(d)
    st = reaction_case(d, False)
    se = reaction_case(d, True)
    st['FS'] = cap['Qamm'] / np.maximum(np.abs(st['R']), 1e-9)
    se['FS'] = cap['Qamm'] / np.maximum(np.abs(se['R']), 1e-9)
    return {
        'stratigrafia': cap['stratigrafia'],
        'Qult_palo': cap['Qult'],
        'Qamm_palo': cap['Qamm'],
        'Qb': cap['Qb'],
        'Qs': cap['Qs'],
        'statico': st,
        'sismico': se,
    }


def tabella_sintesi(r: Dict[str, object]) -> pd.DataFrame:
    rows = [
        ('Qult palo [kN]', r['Qult_palo'], r['Qult_palo']),
        ('Qamm palo [kN]', r['Qamm_palo'], r['Qamm_palo']),
        ('Rmax palo [kN]', float(np.max(r['statico']['R'])), float(np.max(r['sismico']['R']))),
        ('Rmin palo [kN]', float(np.min(r['statico']['R'])), float(np.min(r['sismico']['R']))),
        ('FS minimo [-]', float(np.min(r['statico']['FS'])), float(np.min(r['sismico']['FS']))),
    ]
    return pd.DataFrame(rows, columns=['Parametro', 'Statico', 'Sismico'])


def tabella_pali(r: Dict[str, object], which='statico') -> pd.DataFrame:
    rr = r[which]
    return pd.DataFrame({
        'x [m]': rr['x'],
        'y [m]': rr['y'],
        'Reazione [kN]': rr['R'],
        'FS [-]': rr['FS'],
    })


def tabella_pali_comparativa(r: Dict[str, object]) -> pd.DataFrame:
    """Tabella unica con colonne comparative statico/sismico, ordinata per reazione sismica decrescente."""
    st = r['statico']
    se = r['sismico']
    n = len(st['x'])
    df = pd.DataFrame({
        'Palo': np.arange(1, n + 1, dtype=int),
        'x [m]': np.round(st['x'], 3),
        'y [m]': np.round(st['y'], 3),
        'R_stat [kN]': np.round(st['R'], 1),
        'FS_stat': np.round(st['FS'], 2),
        'R_sis [kN]': np.round(se['R'], 1),
        'FS_sis': np.round(se['FS'], 2),
    })
    df = df.sort_values('R_sis [kN]', ascending=False).reset_index(drop=True)
    return df


def genera_warning(d: DatiPlintoPali, r: dict) -> List[str]:
    """Genera messaggi di avvertimento basati sui risultati del calcolo."""
    warnings = []

    fs_stat_min = float(np.min(r['statico']['FS']))
    fs_sis_min = float(np.min(r['sismico']['FS']))
    r_stat_min = float(np.min(r['statico']['R']))
    r_sis_min = float(np.min(r['sismico']['R']))
    r_stat_max = float(np.max(r['statico']['R']))
    r_sis_max = float(np.max(r['sismico']['R']))
    qamm = float(r['Qamm_palo'])

    if fs_stat_min < 1.5:
        warnings.append(f"⚠ FS minimo statico basso ({fs_stat_min:.2f} < 1.5) sul palo critico")
    if fs_sis_min < 1.2:
        warnings.append(f"⚠ FS minimo sismico molto basso ({fs_sis_min:.2f} < 1.2)")
    if r_stat_min < 0:
        warnings.append(f"⚠ Palo in trazione nel caso statico ({r_stat_min:.0f} kN): verificare collegamento testa palo")
    if r_sis_min < 0:
        warnings.append(f"⚠ Palo in trazione nel caso sismico ({r_sis_min:.0f} kN)")
    if r_stat_max > qamm:
        warnings.append(f"⚠ Reazione massima statica ({r_stat_max:.0f} kN) supera la portata ammissibile ({qamm:.0f} kN)")
    if r_sis_max > qamm * 1.3:
        warnings.append(f"⚠ Reazione massima sismica ({r_sis_max:.0f} kN) supera 1.3×Qamm ({qamm * 1.3:.0f} kN)")
    if d.interasse_x < 3 * d.diametro_palo:
        warnings.append(f"ℹ Interasse x ({d.interasse_x:.2f} m) < 3D ({3 * d.diametro_palo:.2f} m): valutare effetto gruppo")
    if d.interasse_y < 3 * d.diametro_palo:
        warnings.append(f"ℹ Interasse y ({d.interasse_y:.2f} m) < 3D ({3 * d.diametro_palo:.2f} m): valutare effetto gruppo")

    return warnings


def genera_note(d: DatiPlintoPali, r: dict) -> List[str]:
    """Genera note tecniche interpretative sul progetto."""
    note = []
    n_tot = d.n_x * d.n_y
    qamm = float(r['Qamm_palo'])
    qb = float(r.get('Qb', 0.0))
    qs = float(r.get('Qs', 0.0))

    note.append(f"Numero totale pali: {d.n_x}×{d.n_y} = {n_tot}")
    note.append("Reazioni calcolate con distribuzione rigida del plinto (piano rigido)")

    if qs > 0 and qb > 0:
        if qb > qs:
            note.append("La portata di punta governa sulla portata laterale")
        elif qs > qb * 2:
            note.append("La portata laterale di fusto governa sulla punta")

    r_stat = r['statico']['R']
    r_sis = r['sismico']['R']
    idx_max = int(np.argmax(r_stat))
    idx_min = int(np.argmin(r_stat))
    note.append(f"Palo più caricato (statico): n°{idx_max + 1} con {float(r_stat[idx_max]):.0f} kN")
    note.append(f"Palo meno caricato (statico): n°{idx_min + 1} con {float(r_stat[idx_min]):.0f} kN")

    fs_stat_min = float(np.min(r['statico']['FS']))
    fs_sis_min = float(np.min(r['sismico']['FS']))
    if fs_sis_min < fs_stat_min:
        note.append("Il caso sismico governa il progetto del gruppo pali")

    note.append(f"Portata ammissibile palo = {qamm:.0f} kN (FS = {d.gamma_sicurezza:.1f})")

    return note


def figura_geometria(d: DatiPlintoPali, r: Dict[str, object]) -> go.Figure:
    """Pianta del plinto con pali colorati in base alla reazione statica (heatmap RdYlGn)."""
    rr = r['statico']
    R_vals = rr['R']
    x_pali = rr['x']
    y_pali = rr['y']
    n = len(x_pali)

    fig = go.Figure()

    # Rettangolo plinto
    xs = [-d.B / 2, d.B / 2, d.B / 2, -d.B / 2, -d.B / 2]
    ys = [-d.L / 2, -d.L / 2, d.L / 2, d.L / 2, -d.L / 2]
    fig.add_trace(go.Scatter(
        x=xs, y=ys,
        fill='toself',
        mode='lines',
        fillcolor='rgba(173, 216, 230, 0.3)',
        line=dict(color='steelblue', width=2),
        name='Plinto',
    ))

    # Nocciolo del plinto (B/6 x L/6)
    nk_xs = [-d.B / 6, d.B / 6, d.B / 6, -d.B / 6, -d.B / 6]
    nk_ys = [-d.L / 6, -d.L / 6, d.L / 6, d.L / 6, -d.L / 6]
    fig.add_trace(go.Scatter(
        x=nk_xs, y=nk_ys,
        mode='lines',
        line=dict(color='gray', width=1, dash='dot'),
        name='Nocciolo plinto',
    ))

    # Pali come cerchi (scatter) con heatmap RdYlGn invertita (rosso=alto carico, verde=basso)
    R_max = float(np.max(np.abs(R_vals))) if np.max(np.abs(R_vals)) > 0 else 1.0
    raggio_palo = d.diametro_palo / 2.0
    theta = np.linspace(0, 2 * math.pi, 60)
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)

    # Colore normalizzato 0-1 per ogni palo (maggiore reazione = rosso)
    R_norm = (R_vals - float(np.min(R_vals))) / max(float(np.max(R_vals) - np.min(R_vals)), 1e-9)

    for i in range(n):
        cx, cy = float(x_pali[i]), float(y_pali[i])
        xc = (cx + raggio_palo * cos_t).tolist()
        yc = (cy + raggio_palo * sin_t).tolist()
        # Interpolazione colore RdYlGn: rosso alto, verde basso
        t = float(R_norm[i])
        # RdYlGn: 0=verde(0,128,0), 0.5=giallo(255,255,0), 1=rosso(255,0,0)
        if t < 0.5:
            r_col = int(255 * t * 2)
            g_col = int(128 + 127 * t * 2)
            b_col = 0
        else:
            r_col = 255
            g_col = int(255 * (1 - (t - 0.5) * 2))
            b_col = 0
        fill_color = f'rgba({r_col},{g_col},{b_col},0.8)'
        fig.add_trace(go.Scatter(
            x=xc + [xc[0]],
            y=yc + [yc[0]],
            fill='toself',
            fillcolor=fill_color,
            mode='lines',
            line=dict(color='black', width=1),
            showlegend=(i == 0),
            name='Pali' if i == 0 else '',
            hovertext=f'Palo {i + 1}: {float(R_vals[i]):.0f} kN',
            hoverinfo='text',
        ))

    # Annotazioni numero palo e reazione
    for i in range(n):
        cx, cy = float(x_pali[i]), float(y_pali[i])
        fig.add_annotation(
            x=cx, y=cy,
            text=f"<b>{i + 1}</b><br>{float(R_vals[i]):.0f}",
            showarrow=False,
            font=dict(size=9, color='black'),
            xanchor='center',
            yanchor='middle',
        )

    # Scala colore fittizia per la legenda
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(
            size=10,
            color=[float(np.min(R_vals)), float(np.max(R_vals))],
            colorscale='RdYlGn_r',
            showscale=True,
            colorbar=dict(title='R stat [kN]', x=1.02),
        ),
        showlegend=False,
        hoverinfo='none',
    ))

    fig.update_layout(
        title='Geometria del plinto su pali (reazioni statiche)',
        xaxis_title='x [m]',
        yaxis_title='y [m]',
        template='plotly_white',
        legend=dict(x=0, y=1),
    )
    fig.update_yaxes(scaleanchor='x', scaleratio=1)
    return fig


def figura_output(r: Dict[str, object], which='statico') -> go.Figure:
    rr = r[which]
    sizes = np.clip(np.abs(rr['R']) / max(np.max(np.abs(rr['R'])), 1e-9) * 40, 12, 40)
    fig = go.Figure(data=go.Scatter(
        x=rr['x'], y=rr['y'],
        mode='markers+text',
        text=[f'{v:.0f}' for v in rr['R']],
        textposition='top center',
        marker=dict(
            size=sizes,
            color=rr['R'],
            colorscale='RdBu',
            showscale=True,
            colorbar=dict(title='kN'),
        ),
        name=f'Reazioni {which}',
    ))
    fig.update_layout(
        title=f'Reazioni sui pali - {which}',
        xaxis_title='x [m]',
        yaxis_title='y [m]',
        template='plotly_white',
    )
    fig.update_yaxes(scaleanchor='x', scaleratio=1)
    return fig


def figura_comparativa(r: dict) -> go.Figure:
    """Bar chart verticale affiancato: reazioni statiche e sismiche per palo, con linea Qamm."""
    st = r['statico']
    se = r['sismico']
    n = len(st['x'])
    pali_labels = [f"P{i + 1}" for i in range(n)]
    qamm = float(r['Qamm_palo'])

    fig = go.Figure()

    fig.add_trace(go.Bar(
        name='Statico',
        x=pali_labels,
        y=np.round(st['R'], 1),
        marker_color='steelblue',
        text=[f"{v:.0f}" for v in st['R']],
        textposition='outside',
    ))

    fig.add_trace(go.Bar(
        name='Sismico',
        x=pali_labels,
        y=np.round(se['R'], 1),
        marker_color='darkorange',
        text=[f"{v:.0f}" for v in se['R']],
        textposition='outside',
    ))

    # Linea orizzontale Qamm
    fig.add_hline(
        y=qamm,
        line_dash='dash',
        line_color='red',
        annotation_text=f"Qamm = {qamm:.0f} kN",
        annotation_position='top right',
        annotation_font_color='red',
    )

    fig.update_layout(
        title='Confronto reazioni sui pali: caso statico vs sismico',
        xaxis_title='Palo',
        yaxis_title='Reazione [kN]',
        barmode='group',
        template='plotly_white',
        legend=dict(x=0.01, y=0.99),
    )
    return fig
