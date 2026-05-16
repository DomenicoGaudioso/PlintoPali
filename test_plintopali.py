# -*- coding: utf-8 -*-
"""
Test del solver PlintoPali — ripartizione Navier su plinto rigido, capacità palo.

Verifica la ripartizione delle forze su gruppo di pali (ipotesi plinto rigido)
e la capacità portante del singolo palo, senza richiedere OpenSeesPy.
OpenSeesPy viene rimpiazzato da un mock minimale per poter importare src.py.
Dati conformi NTC2018, geometrie realistiche.
"""
import math
import sys
import os
import types

# --- Mock minimale di openseespy per importare src.py senza il modulo installato ---
_ops_mock = types.ModuleType("openseespy")
_ops_sub = types.ModuleType("openseespy.opensees")

def _noop(*args, **kwargs):
    return None

for _fn in [
    "wipe", "model", "node", "fix", "uniaxialMaterial", "element",
    "geomTransf", "equalDOF", "timeSeries", "pattern", "load",
    "system", "numberer", "constraints", "integrator", "algorithm",
    "analysis", "analyze", "basicForce", "nodeDisp", "eleForce",
]:
    setattr(_ops_sub, _fn, _noop)

_ops_mock.opensees = _ops_sub
sys.modules["openseespy"] = _ops_mock
sys.modules["openseespy.opensees"] = _ops_sub
# ---------------------------------------------------------------------------------

sys.path.insert(0, os.path.dirname(__file__))
from src import (
    DatiPlintoPali,
    coordinate_pali,
    reaction_case_rigid,
    capacita_singolo_palo,
    efficienza_gruppo_converse_labarre,
)

import numpy as np


STRAT_SABBIA = "8.0,18,20,30,0,25000"
STRAT_ARGILLA = "8.0,18,20,0,120,60000"


def _dati_4pali(N=1000.0, Mx=0.0, My=0.0, strat=STRAT_SABBIA):
    """Crea DatiPlintoPali con 4 pali simmetrici 2×2."""
    return DatiPlintoPali(
        B=4.0, L=4.0, spessore_plinto=1.0, E_cls_MPa=30000.0,
        n_x=2, n_y=2, interasse_x=2.0, interasse_y=2.0,
        diametro_palo=0.5, lunghezza_palo=8.0,
        Nq=15.0, Nc=9.0, beta=0.3, alpha=0.5,
        gamma_sicurezza=2.5,
        N=N, Mx=Mx, My=My,
        kh=0.1, kv=0.05, falda=100.0,
        stratigrafia_csv=strat,
    )


def test_coordinata_4pali_simmetrici():
    """4 pali in griglia 2×2 con interasse=2m: coordinate centrate in (0,0).

    Le coordinate devono essere simmetriche: ±1 m in x e y.
    """
    d = _dati_4pali()
    x, y = coordinate_pali(d)
    assert len(x) == 4, f"Numero pali = {len(x)}, atteso 4"
    for xi, yi in zip(x, y):
        assert abs(abs(xi) - 1.0) < 0.001, f"x_palo = {xi:.3f}, atteso ±1.0"
        assert abs(abs(yi) - 1.0) < 0.001, f"y_palo = {yi:.3f}, atteso ±1.0"


def test_ripartizione_navier_carico_centrico():
    """N=1000 kN centrico su 4 pali simmetrici → R_i = 250 kN per ogni palo.

    Con e=0 e nessun momento la distribuzione è uniforme per l'ipotesi di plinto rigido.
    Formula Navier: Ri = N/n + Mx×yi/Σyi² + My×xi/Σxi².
    """
    d = _dati_4pali(N=1000.0, Mx=0.0, My=0.0)
    rig = reaction_case_rigid(d, seismic=False)
    R = rig["R"]
    assert len(R) == 4, f"Numero reazioni = {len(R)}, atteso 4"
    for i, r in enumerate(R):
        assert abs(r - 250.0) < 1.0, f"R_palo_{i+1} = {r:.1f} kN, atteso 250.0 kN"


def test_ripartizione_navier_con_momento_mx():
    """N=1000 kN + Mx=200 kN·m su 4 pali simmetrici, interasse_y=2m.

    Σy² = 4 × 1² = 4 m².
    R_max = 250 + 200×1/4 = 300 kN; R_min = 200 kN.
    """
    d = _dati_4pali(N=1000.0, Mx=200.0, My=0.0)
    rig = reaction_case_rigid(d, seismic=False)
    R = rig["R"]
    R_max = float(np.max(R))
    R_min = float(np.min(R))
    assert abs(R_max - 300.0) < 1.0, f"R_max = {R_max:.1f} kN, atteso 300.0 kN"
    assert abs(R_min - 200.0) < 1.0, f"R_min = {R_min:.1f} kN, atteso 200.0 kN"


def test_equilibrio_somma_reazioni():
    """La somma delle reazioni deve essere uguale al carico totale N.

    Equilibrio globale verticale: ΣRi = N (indipendente dalla distribuzione).
    """
    d = _dati_4pali(N=1500.0, Mx=150.0, My=100.0)
    rig = reaction_case_rigid(d, seismic=False)
    R_tot = float(np.sum(rig["R"]))
    N_design = 1500.0 * (1.0 - 0.0)  # kv=0 in caso statico
    assert abs(R_tot - N_design) < 1.0, (
        f"ΣR = {R_tot:.1f} kN ≠ N_design = {N_design:.1f} kN"
    )


def test_capacita_palo_argilla_alpha():
    """Capacità portante palo in argilla cu=120 kPa, metodo α=0.5.

    D=0.5m, L=8m: Qs = 0.5×120×π×0.5×8 ≈ 754 kN.
    Qb = 9×120×π×(0.25)² ≈ 212 kN → Qult ≈ 966 kN.
    """
    d_argilla = DatiPlintoPali(
        B=4.0, L=4.0, spessore_plinto=1.0, E_cls_MPa=30000.0,
        n_x=2, n_y=2, interasse_x=2.0, interasse_y=2.0,
        diametro_palo=0.5, lunghezza_palo=8.0,
        Nq=9.0, Nc=9.0, beta=0.3, alpha=0.5,
        gamma_sicurezza=2.5,
        N=1000.0, Mx=0.0, My=0.0,
        kh=0.1, kv=0.05, falda=100.0,
        stratigrafia_csv=STRAT_ARGILLA,
    )
    cap = capacita_singolo_palo(d_argilla)
    Qs_atteso = 0.5 * 120.0 * math.pi * 0.5 * 8.0
    Qb_atteso = 9.0 * 120.0 * math.pi * (0.25 ** 2)
    assert abs(cap["Qs"] - Qs_atteso) < 2.0, (
        f"Qs = {cap['Qs']:.1f} kN, atteso {Qs_atteso:.1f} kN"
    )
    assert abs(cap["Qb"] - Qb_atteso) < 2.0, (
        f"Qb = {cap['Qb']:.1f} kN, atteso {Qb_atteso:.1f} kN"
    )


def test_effetto_gruppo_pali_stretti():
    """Effetto gruppo (Converse-Labarre): s/D=2.0 < 4 → η < 1.

    Pali troppo vicini riducono la capacità del gruppo.
    """
    d_stretti = DatiPlintoPali(
        B=2.0, L=2.0, spessore_plinto=0.8, E_cls_MPa=30000.0,
        n_x=2, n_y=2, interasse_x=1.0, interasse_y=1.0,
        diametro_palo=0.5, lunghezza_palo=8.0,
        Nq=15.0, Nc=9.0, beta=0.3, alpha=0.5,
        gamma_sicurezza=2.5,
        N=800.0, Mx=0.0, My=0.0,
        kh=0.1, kv=0.05, falda=100.0,
        stratigrafia_csv=STRAT_SABBIA,
    )
    eff = efficienza_gruppo_converse_labarre(d_stretti)
    assert eff["eta"] < 1.0, (
        f"Efficienza = {eff['eta']:.3f} deve essere < 1 per s/D={eff['s_su_D']:.1f}"
    )


def test_effetto_gruppo_pali_larghi():
    """Nessun effetto gruppo per s/D=6.0 > 4 → η = 1.

    Pali ben spaziati lavorano come singoli elementi indipendenti.
    """
    d_larghi = DatiPlintoPali(
        B=8.0, L=8.0, spessore_plinto=1.5, E_cls_MPa=30000.0,
        n_x=2, n_y=2, interasse_x=3.0, interasse_y=3.0,
        diametro_palo=0.5, lunghezza_palo=8.0,
        Nq=15.0, Nc=9.0, beta=0.3, alpha=0.5,
        gamma_sicurezza=2.5,
        N=1000.0, Mx=0.0, My=0.0,
        kh=0.1, kv=0.05, falda=100.0,
        stratigrafia_csv=STRAT_SABBIA,
    )
    eff = efficienza_gruppo_converse_labarre(d_larghi)
    assert abs(eff["eta"] - 1.0) < 0.001, (
        f"Efficienza = {eff['eta']:.3f} deve essere 1.0 per s/D={eff['s_su_D']:.1f}"
    )


if __name__ == "__main__":
    tests = [
        test_coordinata_4pali_simmetrici,
        test_ripartizione_navier_carico_centrico,
        test_ripartizione_navier_con_momento_mx,
        test_equilibrio_somma_reazioni,
        test_capacita_palo_argilla_alpha,
        test_effetto_gruppo_pali_stretti,
        test_effetto_gruppo_pali_larghi,
    ]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"  PASS  {t.__name__}")
        except AssertionError as e:
            print(f"  FAIL  {t.__name__}: {e}")
            failed += 1
        except Exception as e:
            print(f"  ERROR {t.__name__}: {type(e).__name__}: {e}")
            failed += 1
    print(f"\n{len(tests) - failed}/{len(tests)} test superati")
    sys.exit(failed)
