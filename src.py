# -*- coding: utf-8 -*-
from __future__ import annotations
from dataclasses import asdict, dataclass
from typing import Tuple, List, Dict
import json
import math
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import openseespy.opensees as ops

GAMMA_W = 9.81
DEFAULT_STRAT = """2.0,18,20,30,0,25000
3.0,19,21,34,0,40000
5.0,20,20,0,120,60000
"""

# ==========================================
# 1. PARSING E UTILITÀ GEOTECNICHE
# ==========================================

def parse_stratigrafia(csv_text: str) -> Tuple[pd.DataFrame, List[str]]:
    """Righe: spessore,gamma_dry,gamma_sat,phi_deg,cu_kPa,E_ed_kPa"""
    err, rows = [], []
    lines = [ln.strip() for ln in csv_text.splitlines() if ln.strip()]
    
    if not lines:
        return pd.DataFrame(columns=['spessore_m', 'gamma_dry', 'gamma_sat', 'phi_deg', 'cu_kPa', 'E_ed_kPa']), ['Inserire almeno uno strato.']
        
    for i, line in enumerate(lines, start=1):
        parts = [p.strip() for p in line.replace(';', ',').split(',') if p.strip()]
        if len(parts) != 6:
            err.append(f'Riga {i}: usare 6 campi = spessore,gamma_dry,gamma_sat,phi,cu,E_ed.')
            continue
        try:
            h, gd, gs, phi, cu, e_ed = map(float, parts)
            rows.append({
                'spessore_m': h, 
                'gamma_dry': gd, 
                'gamma_sat': gs, 
                'phi_deg': phi, 
                'cu_kPa': cu, 
                'E_ed_kPa': e_ed
            })
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

# ==========================================
# 2. DEFINIZIONE DATI E VALIDAZIONE
# ==========================================

@dataclass(frozen=True)
class DatiPlintoPali:
    B: float
    L: float
    spessore_plinto: float
    E_cls_MPa: float
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
    if d.spessore_plinto <= 0:
        err.append('Lo spessore del plinto deve essere positivo.')
    if d.E_cls_MPa <= 0:
        err.append('Il modulo elastico del calcestruzzo deve essere positivo.')
    if d.n_x <= 0 or d.n_y <= 0:
        err.append('Il numero di pali per direzione deve essere positivo.')
    if d.interasse_x <= 0 or d.interasse_y <= 0:
        err.append('Gli interassi devono essere positivi.')
    if d.lunghezza_palo <= 0 or d.diametro_palo <= 0:
        err.append('Le dimensioni del palo devono essere positive.')
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

# ==========================================
# 3. GEOTECNICA AVANZATA (CAPACITÀ E CEDIMENTI)
# ==========================================

def stima_kv_palo(d: DatiPlintoPali, df: pd.DataFrame) -> float:
    """Stima la rigidezza assiale combinata (strutturale + geotecnica)."""
    E_c = d.E_cls_MPa * 1000 
    A_p = math.pi * (d.diametro_palo**2) / 4.0
    k_strut = (E_c * A_p) / d.lunghezza_palo
    
    sum_E_dz = 0.0
    for _, r in df.iterrows():
        z_t = max(0.0, float(r['z_top_m']))
        z_b = min(d.lunghezza_palo, float(r['z_bot_m']))
        if z_b > z_t:
            sum_E_dz += float(r['E_ed_kPa']) * (z_b - z_t)
            
    E_s_avg = sum_E_dz / d.lunghezza_palo if d.lunghezza_palo > 0 else 10000.0
    k_geo = E_s_avg * d.diametro_palo * 2.0 
    
    if k_geo <= 0: return k_strut
    return 1.0 / (1.0 / k_strut + 1.0 / k_geo)

def efficienza_gruppo_converse_labarre(d: DatiPlintoPali) -> dict:
    """Valuta l'effetto gruppo secondo Converse-Labarre."""
    s_medio = (d.interasse_x + d.interasse_y) / 2.0
    rapporto_s_D = s_medio / d.diametro_palo
    
    if rapporto_s_D > 4.0 or (d.n_x == 1 and d.n_y == 1):
        return {'eta': 1.0, 'stato': 'Pali Singoli', 's_su_D': rapporto_s_D}
        
    theta = math.degrees(math.atan(d.diametro_palo / s_medio))
    m, n = d.n_x, d.n_y
    eta = 1.0 - (theta / 90.0) * (((n - 1) * m + (m - 1) * n) / (m * n))
    eta = max(0.1, min(1.0, eta))
    
    stato = "Comportamento a Gruppo" if rapporto_s_D <= 3.0 else "Interferenza Parziale"
    return {'eta': eta, 'stato': stato, 's_su_D': rapporto_s_D}

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
        zmid = 0.5 * (a + b)
        phi, cu = float(r['phi_deg']), float(r['cu_kPa'])
        qs = d.beta * sigma_v_eff(df, zmid, d.falda) if phi > 0 else d.alpha * cu
        Qs += qs * P * (b - a)
        
    tip = layer_at_depth(df, d.lunghezza_palo)
    if float(tip['phi_deg']) > 0:
        Qb = d.Nq * sigma_v_eff(df, d.lunghezza_palo, d.falda) * A_base
    else:
        Qb = d.Nc * float(tip['cu_kPa']) * A_base
        
    Qult = Qb + Qs
    return {'stratigrafia': df, 'Qult': Qult, 'Qamm': Qult / d.gamma_sicurezza, 'Qb': Qb, 'Qs': Qs}

def calcola_cedimento_gruppo(d: DatiPlintoPali, Q_tot: float) -> float:
    """Metodo della Zattera Equivalente a 2/3 L."""
    df, _ = parse_stratigrafia(d.stratigrafia_csv)
    z_eq = (2.0 / 3.0) * d.lunghezza_palo
    B_eq = (d.n_x - 1) * d.interasse_x + d.diametro_palo
    L_eq = (d.n_y - 1) * d.interasse_y + d.diametro_palo
    cedimento_m = 0.0
    
    for _, r in df.iterrows():
        z_top = max(z_eq, float(r['z_top_m']))
        z_bot = max(z_eq, float(r['z_bot_m']))
        if z_bot <= z_top: continue
        dz = z_bot - z_top
        z_mid = z_top + dz/2.0 - z_eq
        
        q_z = Q_tot / ((B_eq + z_mid) * (L_eq + z_mid))
        if float(r['E_ed_kPa']) > 0:
            cedimento_m += (q_z / float(r['E_ed_kPa'])) * dz
            
    return cedimento_m * 1000.0

# ==========================================
# 4. MOTORI STRUTTURALI: RIGIDO E FEM
# ==========================================

def reaction_case_rigid(d: DatiPlintoPali, seismic=False):
    """Calcolo reazioni con ipotesi di plinto infinitamente rigido (Navier)."""
    x, y = coordinate_pali(d)
    n = len(x)
    A = np.vstack([np.ones(n), y, x]).T
    N = d.N * (1.0 - (d.kv if seismic else 0.0))
    Mx = d.Mx * (1.0 + (d.kh if seismic else 0.0))
    My = d.My * (1.0 + (d.kh if seismic else 0.0))
    
    abc = np.linalg.solve(A.T @ A, np.array([N, Mx, My]))
    R = A @ abc
    return {'x': x, 'y': y, 'R': R}

def calcola_strut_and_tie(d: DatiPlintoPali, rig_res: dict) -> dict:
    """Calcolo tiranti per plinto tozzo."""
    d_utile = d.spessore_plinto - 0.10
    if d_utile <= 0: d_utile = 0.1
    x, y, R = rig_res['x'], rig_res['y'], rig_res['R']
    
    Mx_tie = sum(R[i] * abs(x[i]) for i in range(len(R)) if R[i] > 0)
    My_tie = sum(R[i] * abs(y[i]) for i in range(len(R)) if R[i] > 0)
    
    Tx = Mx_tie / (0.9 * d_utile)
    Ty = My_tie / (0.9 * d_utile)
    return {'Tx_kN': Tx, 'Ty_kN': Ty, 'd_utile_m': d_utile}

def reaction_case_fem(d: DatiPlintoPali, k_v_calcolato: float, seismic=False):
    """Calcolo reazioni e deformazioni tramite OpenSeesPy con check lunghezza nulla."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    x, y = coordinate_pali(d)
    n = len(x)
    
    E_c = d.E_cls_MPa * 1e3
    G_c = E_c / (2 * (1 + 0.2))
    A_beam = d.interasse_x * d.spessore_plinto 
    I_beam = (d.interasse_x * d.spessore_plinto**3) / 12.0
    J_tors = I_beam * 2 
    
    ops.geomTransf('Linear', 1, 0, 0, 1)
    
    # Nodo Baricentro (Tag 1000)
    ops.node(1000, 0.0, 0.0, 0.0) 
    
    for i in range(n):
        node_tag = i + 1
        px, py = float(x[i]), float(y[i])
        
        # Calcolo distanza dal baricentro
        dist = math.sqrt(px**2 + py**2)
        
        # Nodo Palo (Tag 1..n)
        ops.node(node_tag, px, py, 0.0)
        
        # Nodo Base Molla (Tag 2000..n)
        ops.node(2000+i, px, py, 0.0)
        ops.fix(2000+i, 1, 1, 1, 1, 1, 1)
        
        # Materiale e Molla
        ops.uniaxialMaterial('Elastic', node_tag, k_v_calcolato)
        ops.element('zeroLength', node_tag, 2000+i, node_tag, '-mat', node_tag, '-dir', 3)
        
        # --- FIX: Controllo lunghezza elemento ---
        if dist < 1e-4:
            # Se il palo è nel baricentro, vincoliamo i nodi insieme (Rigid Link)
            ops.equalDOF(1000, node_tag, 1, 2, 3, 4, 5, 6)
        else:
            # Altrimenti creiamo la trave del graticcio
            ops.element('elasticBeamColumn', 1000+i, 1000, node_tag, A_beam, E_c, G_c, J_tors, I_beam, I_beam, 1)

    # Vincoli globali baricentro
    ops.fix(1000, 1, 1, 0, 0, 0, 1)
    
    # Carichi
    N_load = -d.N * (1.0 - (d.kv if seismic else 0.0))
    Mx_load = d.Mx * (1.0 + (d.kh if seismic else 0.0))
    My_load = d.My * (1.0 + (d.kh if seismic else 0.0))
    
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(1000, 0.0, 0.0, N_load, Mx_load, My_load, 0.0)
    
    ops.system('BandSPD')
    ops.numberer('RCM')
    ops.constraints('Transformation') # Cambiato in Transformation per gestire equalDOF
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')
    
    if ops.analyze(1) != 0:
        raise RuntimeError("OpenSees non è riuscito a convergere. Verificare labilità del modello.")
    
    R_fem, settlements, M_radice = np.zeros(n), np.zeros(n), np.zeros(n)
    for i in range(n):
        node_tag = i + 1
        R_fem[i] = -ops.basicForce(node_tag)[0]
        settlements[i] = -ops.nodeDisp(node_tag, 3) * 1000.0
        
        # Recupero momento solo se l'elemento esiste (dist > 0)
        dist = math.sqrt(float(x[i])**2 + float(y[i])**2)
        if dist >= 1e-4:
            ele_forces = ops.eleForce(1000+i)
            M_radice[i] = math.sqrt(ele_forces[4]**2 + ele_forces[5]**2)
        else:
            M_radice[i] = 0.0
            
    return {'x': x, 'y': y, 'R': R_fem, 'cedimenti_mm': settlements, 'M_radice_kNm': M_radice}


def calcola_plinto_pali(d: DatiPlintoPali) -> Dict[str, object]:
    """Orchestratore principale dell'analisi."""
    cap = capacita_singolo_palo(d)
    k_v_calc = stima_kv_palo(d, cap['stratigrafia'])
    eff_gruppo = efficienza_gruppo_converse_labarre(d)
    Qamm_effettiva = cap['Qamm'] * eff_gruppo['eta']
    
    st_rig = reaction_case_rigid(d, False)
    se_rig = reaction_case_rigid(d, True)
    
    st_fem = reaction_case_fem(d, k_v_calc, False)
    se_fem = reaction_case_fem(d, k_v_calc, True)
    
    st_rig['FS'] = cap['Qamm'] / np.maximum(np.abs(st_rig['R']), 1e-9)
    se_rig['FS'] = cap['Qamm'] / np.maximum(np.abs(se_rig['R']), 1e-9)
    st_fem['FS'] = cap['Qamm'] / np.maximum(np.abs(st_fem['R']), 1e-9)
    se_fem['FS'] = cap['Qamm'] / np.maximum(np.abs(se_fem['R']), 1e-9)
    
    tie_forces = calcola_strut_and_tie(d, st_rig)
    ced_gruppo = calcola_cedimento_gruppo(d, d.N)
    
    return {
        'stratigrafia': cap['stratigrafia'],
        'Qult_palo': cap['Qult'],
        'Qamm_palo': cap['Qamm'],
        'Qb': cap['Qb'],
        'Qs': cap['Qs'],
        'statico': st_rig,
        'sismico': se_rig,
        'statico_fem': st_fem,
        'sismico_fem': se_fem,
        'k_v_calcolato': k_v_calc,
        'efficienza_gruppo': eff_gruppo,
        'Qamm_effettiva_palo': Qamm_effettiva,
        'cedimento_gruppo_mm': ced_gruppo,
        'strut_and_tie': tie_forces
    }

# ==========================================
# 5. REPORTISTICA TABELLARE E WARNING
# ==========================================

def tabella_sintesi(r: Dict[str, object]) -> pd.DataFrame:
    rows = [
        ('Qult palo isolato [kN]', r['Qult_palo'], r['Qult_palo']),
        ('Qamm (effetto gruppo) [kN]', r['Qamm_effettiva_palo'], r['Qamm_effettiva_palo']),
        ('Rmax palo (Rigido) [kN]', float(np.max(r['statico']['R'])), float(np.max(r['sismico']['R']))),
        ('Rmin palo (Rigido) [kN]', float(np.min(r['statico']['R'])), float(np.min(r['sismico']['R']))),
        ('Rmax palo (FEM) [kN]', float(np.max(r['statico_fem']['R'])), float(np.max(r['sismico_fem']['R']))),
        ('FS minimo (Rigido) [-]', float(np.min(r['statico']['FS'])), float(np.min(r['sismico']['FS']))),
    ]
    return pd.DataFrame(rows, columns=['Parametro', 'Statico', 'Sismico'])

def tabella_pali_comparativa(r: Dict[str, object]) -> pd.DataFrame:
    st = r['statico']
    se = r['sismico']
    st_fem = r['statico_fem']
    n = len(st['x'])
    
    df = pd.DataFrame({
        'Palo': np.arange(1, n + 1, dtype=int),
        'x [m]': np.round(st['x'], 3),
        'y [m]': np.round(st['y'], 3),
        'R_Rigido_Stat [kN]': np.round(st['R'], 1),
        'R_FEM_Stat [kN]': np.round(st_fem['R'], 1),
        'Cedimento_FEM [mm]': np.round(st_fem['cedimenti_mm'], 2),
        'R_Rigido_Sis [kN]': np.round(se['R'], 1),
        'FS_Rigido_Stat': np.round(st['FS'], 2),
    })
    return df.sort_values('R_Rigido_Stat [kN]', ascending=False).reset_index(drop=True)

def genera_warning(d: DatiPlintoPali, r: dict) -> List[str]:
    warnings = []
    fs_stat_min = float(np.min(r['statico']['FS']))
    fs_sis_min = float(np.min(r['sismico']['FS']))
    r_stat_min = float(np.min(r['statico']['R']))
    r_sis_min = float(np.min(r['sismico']['R']))
    r_stat_max = float(np.max(r['statico']['R']))
    r_fem_max = float(np.max(r['statico_fem']['R']))
    qamm_eff = float(r['Qamm_effettiva_palo'])

    if fs_stat_min < 1.5:
        warnings.append(f"⚠ FS minimo statico basso ({fs_stat_min:.2f} < 1.5)")
    if fs_sis_min < 1.2:
        warnings.append(f"⚠ FS minimo sismico molto basso ({fs_sis_min:.2f} < 1.2)")
    if r_stat_min < 0:
        warnings.append(f"⚠ Palo in trazione nel caso statico ({r_stat_min:.0f} kN)")
    if r_sis_min < 0:
        warnings.append(f"⚠ Palo in trazione nel caso sismico ({r_sis_min:.0f} kN)")
    if r_stat_max > qamm_eff:
        warnings.append(f"⚠ Reazione massima rigida ({r_stat_max:.0f} kN) supera Qamm efficace ({qamm_eff:.0f} kN)")
    if r_fem_max > qamm_eff:
        warnings.append(f"⚠ Reazione massima FEM ({r_fem_max:.0f} kN) supera Qamm efficace ({qamm_eff:.0f} kN)")
        
    eff = r['efficienza_gruppo']
    if eff['eta'] < 1.0:
        warnings.append(f"ℹ Effetto gruppo attivo (s/D={eff['s_su_D']:.1f}): efficienza ridotta a {eff['eta']:.2f}")

    if d.spessore_plinto < max(d.interasse_x, d.interasse_y) / 2:
        warnings.append("⚠ Plinto sottile rispetto all'interasse: l'ipotesi di plinto rigido potrebbe non essere conservativa.")

    return warnings

def genera_note(d: DatiPlintoPali, r: dict) -> List[str]:
    note = []
    n_tot = d.n_x * d.n_y
    qb, qs = float(r.get('Qb', 0.0)), float(r.get('Qs', 0.0))

    note.append(f"Numero totale pali: {d.n_x}×{d.n_y} = {n_tot}")
    if qs > 0 and qb > 0:
        if qb > qs:
            note.append("La portata di punta governa sulla portata laterale")
        elif qs > qb * 2:
            note.append("La portata laterale di fusto governa sulla punta")

    r_stat = r['statico']['R']
    idx_max = int(np.argmax(r_stat))
    note.append(f"Palo più caricato (ipotesi rigida): n°{idx_max + 1} con {float(r_stat[idx_max]):.0f} kN")
    note.append(f"Modulo assiale k_v calcolato per molle FEM: {r['k_v_calcolato']:.0f} kN/m")

    return note

# ==========================================
# 6. VISUALIZZAZIONI PLOTLY
# ==========================================

def figura_3d_plinto_pali(d: DatiPlintoPali) -> go.Figure:
    fig = go.Figure()
    
    # Plinto
    x0, x1 = -d.B/2, d.B/2
    y0, y1 = -d.L/2, d.L/2
    z0, z1 = 0, d.spessore_plinto
    
    fig.add_trace(go.Mesh3d(
        x=[x0, x0, x1, x1, x0, x0, x1, x1],
        y=[y0, y1, y1, y0, y0, y1, y1, y0],
        z=[z0, z0, z0, z0, z1, z1, z1, z1],
        i=[7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2],
        j=[3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3],
        k=[0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6],
        opacity=0.4, 
        color='steelblue', 
        name='Plinto'
    ))
    
    # Pali
    x_p, y_p = coordinate_pali(d)
    for i in range(len(x_p)):
        fig.add_trace(go.Scatter3d(
            x=[x_p[i], x_p[i]], 
            y=[y_p[i], y_p[i]], 
            z=[0, -d.lunghezza_palo],
            mode='lines', 
            line=dict(color='gray', width=8), 
            name='Pali' if i == 0 else '', 
            showlegend=(i == 0)
        ))
        
    fig.update_layout(
        scene=dict(aspectmode='data', zaxis_title='Z [m]', xaxis_title='X [m]', yaxis_title='Y [m]'), 
        title="Modello 3D Geometria", 
        margin=dict(l=0, r=0, b=0, t=40), 
        template='plotly_white'
    )
    return fig

def figura_stratigrafia(r: dict) -> go.Figure:
    df = r['stratigrafia']
    fig = go.Figure()
    colors = ['#D2B48C', '#DEB887', '#F4A460', '#CD853F', '#A0522D', '#8B4513']
    
    for i, row in df.iterrows():
        z_t, z_b = -row['z_top_m'], -row['z_bot_m']
        c = colors[i % len(colors)]
        fig.add_shape(type="rect", x0=0, y0=z_t, x1=1, y1=z_b, fillcolor=c, line=dict(color="black"))
        testo = f"E_ed={row['E_ed_kPa']:.0f} kPa<br>φ={row['phi_deg']}° cu={row['cu_kPa']} kPa"
        fig.add_annotation(x=0.5, y=(z_t+z_b)/2, text=testo, showarrow=False, font=dict(color="black", size=10))
        
    fig.update_layout(
        title="Profilo Stratigrafico", 
        xaxis=dict(visible=False), 
        yaxis_title="Profondità Z [m]", 
        template='plotly_white'
    )
    return fig

def figura_mesh_fem(d: DatiPlintoPali, r: dict) -> go.Figure:
    fig = go.Figure()
    st_fem = r['statico_fem']
    x_p, y_p, M_rad = st_fem['x'], st_fem['y'], st_fem['M_radice_kNm']
    
    # Nodi alla base e colore del cedimento
    fig.add_trace(go.Scatter3d(
        x=x_p, y=y_p, z=[0]*len(x_p), 
        mode='markers+text',
        marker=dict(
            size=6, 
            color=st_fem['cedimenti_mm'], 
            colorscale='Viridis', 
            showscale=True, 
            colorbar=dict(title="Cedimento [mm]")
        ),
        text=[f"P{i+1}" for i in range(len(x_p))], 
        textposition="bottom center", 
        name='Nodi Pali'
    ))
    
    # Baricentro
    fig.add_trace(go.Scatter3d(
        x=[0], y=[0], z=[0], 
        mode='markers', 
        marker=dict(size=8, color='red', symbol='cross'), 
        name='Baricentro'
    ))
    
    # Graticcio
    for i in range(len(x_p)):
        fig.add_trace(go.Scatter3d(
            x=[0, x_p[i]], y=[0, y_p[i]], z=[0, 0], 
            mode='lines',
            line=dict(color='darkred', width=4), 
            hovertext=f"Momento Radice: {M_rad[i]:.1f} kNm", 
            hoverinfo="text", 
            name='Graticcio FEM' if i == 0 else '',
            showlegend=(i == 0)
        ))
        
    fig.update_layout(
        scene=dict(zaxis=dict(range=[-1, 1], visible=False)), 
        title="Mesh FEM Graticcio", 
        margin=dict(l=0, r=0, b=0, t=40), 
        template='plotly_white'
    )
    return fig

def figura_geometria(d: DatiPlintoPali, r: Dict[str, object]) -> go.Figure:
    rr = r['statico']
    R_vals, x_pali, y_pali = rr['R'], rr['x'], rr['y']
    fig = go.Figure()
    
    xs = [-d.B / 2, d.B / 2, d.B / 2, -d.B / 2, -d.B / 2]
    ys = [-d.L / 2, -d.L / 2, d.L / 2, d.L / 2, -d.L / 2]
    
    fig.add_trace(go.Scatter(
        x=xs, y=ys, 
        fill='toself', 
        mode='lines', 
        fillcolor='rgba(173, 216, 230, 0.3)', 
        line=dict(color='steelblue', width=2), 
        name='Plinto'
    ))

    R_max_val = np.max(R_vals) if np.max(R_vals) > np.min(R_vals) else np.max(R_vals) + 1.0
    R_norm = (R_vals - float(np.min(R_vals))) / float(R_max_val - np.min(R_vals))
    theta = np.linspace(0, 2 * math.pi, 30)
    
    for i in range(len(x_pali)):
        xc = (x_pali[i] + (d.diametro_palo/2) * np.cos(theta)).tolist()
        yc = (y_pali[i] + (d.diametro_palo/2) * np.sin(theta)).tolist()
        t = float(R_norm[i])
        
        if t < 0.5:
            r_col, g_col, b_col = int(255*t*2), int(128+127*t*2), 0
        else:
            r_col, g_col, b_col = 255, int(255*(1-(t-0.5)*2)), 0
            
        fig.add_trace(go.Scatter(
            x=xc+[xc[0]], y=yc+[yc[0]], 
            fill='toself', 
            fillcolor=f'rgba({r_col},{g_col},{b_col},0.8)', 
            mode='lines', 
            line=dict(color='black', width=1), 
            showlegend=False, 
            hovertext=f'Palo {i+1}: {float(R_vals[i]):.0f} kN'
        ))
        
        fig.add_annotation(
            x=x_pali[i], y=y_pali[i], 
            text=f"<b>{i+1}</b><br>{float(R_vals[i]):.0f}", 
            showarrow=False, 
            font=dict(size=9)
        )
        
    fig.update_layout(
        title='Geometria in Pianta (Reazioni Statiche Rigide)', 
        template='plotly_white', 
        yaxis=dict(scaleanchor='x', scaleratio=1)
    )
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
    st = r['statico']
    se = r['sismico']
    pali_labels = [f"P{i + 1}" for i in range(len(st['x']))]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        name='Statico', x=pali_labels, y=np.round(st['R'], 1), 
        marker_color='steelblue', text=[f"{v:.0f}" for v in st['R']], textposition='outside'
    ))
    fig.add_trace(go.Bar(
        name='Sismico', x=pali_labels, y=np.round(se['R'], 1), 
        marker_color='darkorange', text=[f"{v:.0f}" for v in se['R']], textposition='outside'
    ))

    qamm = float(r['Qamm_effettiva_palo'])
    fig.add_hline(
        y=qamm, line_dash='dash', line_color='red', 
        annotation_text=f"Qamm Eff. = {qamm:.0f} kN"
    )

    fig.update_layout(
        title='Confronto reazioni: Statico vs Sismico (Ipotesi Rigida)',
        barmode='group', 
        template='plotly_white'
    )
    return fig

def figura_comparativa_rigido_fem(r: dict) -> go.Figure:
    st_rig, st_fem = r['statico'], r['statico_fem']
    pali_labels = [f"P{i + 1}" for i in range(len(st_rig['x']))]
    
    fig = go.Figure()
    fig.add_trace(go.Bar(
        name='Plinto Rigido', x=pali_labels, y=np.round(st_rig['R'], 1), 
        marker_color='steelblue', text=[f"{v:.0f}" for v in st_rig['R']], textposition='outside'
    ))
    fig.add_trace(go.Bar(
        name='Plinto Flessibile (FEM)', x=pali_labels, y=np.round(st_fem['R'], 1), 
        marker_color='darkorange', text=[f"{v:.0f}" for v in st_fem['R']], textposition='outside'
    ))
    
    qamm = float(r['Qamm_effettiva_palo'])
    fig.add_hline(
        y=qamm, line_dash='dash', line_color='red', 
        annotation_text=f"Qamm Eff. = {qamm:.0f} kN"
    )
    
    fig.update_layout(
        title='Confronto Reazioni Statiche: Plinto Rigido vs Flessibile', 
        barmode='group', 
        template='plotly_white'
    )
    return fig