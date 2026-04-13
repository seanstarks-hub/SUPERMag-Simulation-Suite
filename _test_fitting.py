import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
from pathlib import Path
import supermag
# Force Python fallback to test quickly
import supermag.depairing as _dp
_dp._USE_NATIVE = False

nb = supermag.get_material('Nb')
cuni = supermag.get_material('Cu0.43Ni0.57')
Tc0 = 7.1
d_S = 25.0

csv_path = Path('validation/buzdin_1982/expected/nb_cuni_buzdin_figure_9.csv')
raw = np.genfromtxt(csv_path, delimiter=',', skip_header=1)
d_F_data = np.ascontiguousarray(raw[:, 0])
Tc_data  = np.ascontiguousarray(raw[:, 1])
mask = d_F_data > 0
d_F_data = d_F_data[mask]
Tc_data  = Tc_data[mask]
print(f'Data: {len(d_F_data)} pts')

result = supermag.fit_tc(
    Tc0=Tc0, d_S=d_S,
    xi_S=nb['xi_S'], xi_F=cuni['xi_F'], E_ex=cuni['E_ex'],
    gamma=0.1, gamma_B=0.3,
    d_F_data=d_F_data, Tc_data=Tc_data,
    fit_gamma=True, fit_gamma_B=False,
    model='fominov',
)
print(f'gamma = {result["gamma"]:.4f}')
print(f'chi2 = {result["chi2"]:.4f}')

result_2p = supermag.fit_tc(
    Tc0=Tc0, d_S=d_S,
    xi_S=nb['xi_S'], xi_F=cuni['xi_F'], E_ex=cuni['E_ex'],
    gamma=0.1, gamma_B=0.1,
    d_F_data=d_F_data, Tc_data=Tc_data,
    fit_gamma=True, fit_gamma_B=True,
    model='fominov',
)
print(f'2p gamma = {result_2p["gamma"]:.4f}, gamma_B = {result_2p["gamma_B"]:.4f}, chi2 = {result_2p["chi2"]:.4f}')

d_opt = supermag.optimize_tc(
    Tc0=Tc0, d_S=d_S,
    xi_S=nb['xi_S'], xi_F=cuni['xi_F'], E_ex=cuni['E_ex'],
    gamma=0.15, gamma_B=0.3,
    d_F_lo=0.5, d_F_hi=20.0,
    Tc_target=5.0, model='fominov',
)
print(f'optimize_tc d_F = {d_opt:.4f}')

d_exact = supermag.inverse_tc(
    Tc0=Tc0, d_S=d_S,
    xi_S=nb['xi_S'], xi_F=cuni['xi_F'], E_ex=cuni['E_ex'],
    gamma=0.15, gamma_B=0.3,
    Tc_target=5.0, d_F_lo=0.5, d_F_hi=20.0,
    model='fominov',
)
print(f'inverse_tc d_F = {d_exact:.6f}')
print('ALL OK')
