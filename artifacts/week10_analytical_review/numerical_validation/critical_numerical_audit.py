"""Read-only independent 70-digit check of Week 10 population benchmarks."""
from decimal import Decimal as D, localcontext
import hashlib
import json
from pathlib import Path
import numpy as np
from ou_irregular.week8.design import E4Cell
from ou_irregular.week10.theory import population_predictions, laplace_kernel
from ou_irregular.week10.cli import verify_seal


def mul(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def inv(a):
    det = a[0][0]*a[1][1]-a[0][1]*a[1][0]
    return [[a[1][1]/det, -a[0][1]/det], [-a[1][0]/det, a[0][0]/det]]


def independent(ratio, n):
    # Derive directly from model calibration, without project phase_type/kernel.
    with localcontext() as ctx:
        ctx.prec = 70
        r, m = D(ratio), D('.5')
        f = D('.0000005') if ratio == 20 else D(0)
        l0, l1 = D(6)/(r+5)/m, D(6)*r/(r+5)/m
        a, b = D(1)/(10*m), D(1)/(2*m)
        A, L = [[l0+a,-a],[-b,l1+b]], [[l0,D(0)],[D(0),l1]]
        I = [[D(1),D(0)],[D(0),D(1)]]
        tr = (A[0][0]+A[1][1])/2
        d = (((A[0][0]-A[1][1])/2)**2+A[0][1]*A[1][0]).sqrt()
        ep, em = (d*f).exp(), (-d*f).exp()
        ch, sh = (ep+em)/2, (ep-em)/(2*d)
        E = [[(-tr*f).exp()*(ch*I[i][j]-sh*(A[i][j]-tr*I[i][j]))
              for j in range(2)] for i in range(2)]
        head = mul(mul(inv(A), [[I[i][j]-E[i][j] for j in range(2)] for i in range(2)]), L)
        tail = mul(E, mul(inv([[A[i][j]+I[i][j] for j in range(2)] for i in range(2)]), L))
        B = [[(-f).exp()*(head[i][j]+tail[i][j]) for j in range(2)] for i in range(2)]
        v, c = [D(5)/(r+5),r/(r+5)], [D(1)]
        for _ in range(n-1):
            v = [v[0]*B[0][0]+v[1]*B[1][0],v[0]*B[0][1]+v[1]*B[1][1]]
            c.append(sum(v))
        rho, k = c[1], D(n-1)
        variance = (D(n-2)/k**2+(1+rho*rho-2*rho*c[-1])/(k*k*(1-rho)**2)
                    +2*sum(D(n-1-h)*c[h] for h in range(1,n-1))/k**2)/8
        actual = population_predictions(E4Cell(n=1000,q=.5,mechanism='mmpp',ratio=ratio),n)['fixed_rho_mean_variance']
        return {'ratio': ratio, 'n': n, 'decimal_variance': str(variance),
                'implemented_variance': actual, 'relative_error': float((D(actual)-variance)/variance)}

checks = [independent(r,n) for r in [1,5,20] for n in [2,7,1000]]
spectral = []
for r in [1,5,20]:
    kernel = laplace_kernel(E4Cell(n=1000,q=.5,mechanism='mmpp',ratio=r),1.)
    spectral.append({'ratio': r,'spectral_radius': float(max(abs(np.linalg.eigvals(kernel))))})
small = population_predictions(E4Cell(n=1000,q=1e-12,mechanism='iid',ratio=1),n=2)
result = {'status': 'passed', 'scope': 'adopted q=.5 grid; separate unsupported near-unit stress example',
          'arithmetic_precision_decimal_digits': 70,'finite_count_checks':checks,
          'max_relative_error':max(abs(x['relative_error']) for x in checks),
          'kernel_spectral_radius':spectral,
          'off_grid_roundoff_example':{'q':1e-12,'n':2,'mechanism':'iid','ratio':1,
            'relative_variance_error_vs_exact_identity':small['fixed_rho_mean_variance']/small['population_working_variance']-1},
          'current_theory_sha256':hashlib.sha256(Path('ou_irregular/week10/theory.py').read_bytes()).hexdigest(),
          'sealed_package':verify_seal('artifacts/week10_supplement')}
assert result['max_relative_error'] < 1e-14
print(json.dumps(result,indent=2,allow_nan=False))
