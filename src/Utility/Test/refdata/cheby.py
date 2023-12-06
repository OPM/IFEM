import numpy as np

def write(data, name, dim, file):
    file.write(f'const std::string cheby_{name} =\n')
    file.write('R"({} {} {}'.format(dim[0], dim[1], dim[2]))
    for p in data:
        file.write(f'\n{str(p)}')
    file.write(')";\n\n')

xc = [0.0, 0.5, 0.6, 0.8]
u = lambda xi: (xi + 1.0) / 2.0
f1 = lambda x: x**3 + x*x + x - 1.0
f2 = lambda x: x**3 - 2*x*x + 3*x - 1.0
f3 = lambda x: 4*x**3 - 7*x*x - 3*x + 5.0
y1 = [f1(u(x)) for x in xc]
y2 = [f2(u(x)) for x in xc]
y3 = [f3(u(x)) for x in xc]

p1 = np.polynomial.chebyshev.chebfit(xc, y1, deg=3)
q1 = np.polynomial.chebyshev.chebfit(xc, y2, deg=3)
r1 = np.polynomial.chebyshev.chebfit(xc, y3, deg=3)

with open('cheby_strings.h','w') as file:
    write(p1, '1', [4,0,0], file)
    write(np.kron(p1, p1), '2', [4,4,0], file)
    write(np.kron(np.kron(p1, p1), p1), '3', [4,4,4], file)
    write(np.kron(q1, p1), '12', [4,4,0], file)
    write(np.kron(p1, q1), '21', [4,4,0], file)
    write(np.kron(r1, p1), '13', [4,4,0], file)
    write(np.kron(r1, q1), '23', [4,4,0], file)
    write(np.kron(r1, np.kron(q1, p1)), '123', [4,4,4], file)
    write(np.kron(r1, np.kron(p1, q1)), '213', [4,4,4], file)
    write(np.kron(q1, np.kron(r1, p1)), '132', [4,4,4], file)
    write(np.kron(p1, np.kron(r1, q1)), '231', [4,4,4], file)
    write(np.kron(p1, np.kron(q1, r1)), '321', [4,4,4], file)
    write(np.kron(q1, np.kron(p1, r1)), '312', [4,4,4], file)
