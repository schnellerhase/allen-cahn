from fenics import *

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("initial", type=str, choices=['bump', 'dumbel'])
parser.add_argument("--eps", default=.1, type=float)
parser.add_argument("--dimension", default=2, choices=[1, 2], type=int)

args = parser.parse_args()
eps = args.eps
initial_data = args.initial
dimension = args.dimension

dt = 5.0e-5
t = 0.00
T = 0.05

def distance(a, b):
    return sqrt(sum(((a[i] - b[i])**2) for i in range(dimension)))


class LinearBump(UserExpression):
    def eval(self, values, x):
        d = distance(x, (.5, .5))

        if d + DOLFIN_EPS < 0.25:
            values[0] = 1
        elif d + DOLFIN_EPS < 0.25 + eps:
            values[0] = -2 / eps * (d - 0.25 - eps) - 1
        else:
            values[0] = -1

# TODO: fix 1D case
class Dumbel(UserExpression):
    def eval(self, values, x):
        r = 3 / 16
        a = 3 / 8

        if x[0] < 0.5 + DOLFIN_EPS:
            d = distance(x, (a, 0.5))
            if d < r + DOLFIN_EPS:
                values[0] = 1
            elif d < r + eps + DOLFIN_EPS:
                values[0] = 2 / eps * (r - d) + 1
            else:
                values[0] = -1
        else:
            d = distance(x, (1 - a, 0.5))
            if d < r + DOLFIN_EPS:
                values[0] = 1
            elif d < r + eps + DOLFIN_EPS:
                values[0] = 2 / eps * (r - d) + 1
            else:
                values[0] = -1


class PeriodicBoundary(SubDomain):
    def inside(self, x, on_boundary):
        is_near = (near(x[i], 0) for i in range(dimension))
        return bool(any(is_near) and on_boundary)

    def map(self, x, y):
        for i in range(dimension):
            y[i] = x[i] + (near(x[i], 1)) * (-1)


mesh = UnitIntervalMesh(100) if dimension == 1 else UnitSquareMesh(100, 100)

pbc = PeriodicBoundary()
V = FunctionSpace(mesh, "CG", 1, constrained_domain=pbc)

u, v = Function(V), TestFunction(V)
u.rename("u", "")

u_init = LinearBump() if initial_data == "bump" else Dumbel()
u_init = interpolate(u_init, V)

u_pre = Function(V)
u_pre.rename("u", "")
u_pre.interpolate(u_init)


def W(u):
    return (u**2 - 1) ** 2


def W_prime(u):
    return 4 * u * (u**2 - 1)


F = (
    +1 / dt * u * v * dx
    - 1 / dt * u_pre * v * dx
    + inner(grad(u), grad(v)) * dx
    + 1 / eps**2 * W_prime(u) * v * dx
)

i = 0
file = File(f"data/{initial_data}/eps_{eps}.pvd")
file << (u_pre, i)

while t < T:
    print(f"Time {t} of {T}.")

    i += 1
    t += dt
    solve(F == 0, u)
    u_pre.assign(u)

    if i % 10 == 0:
        file << (u_pre, i)


file << (u_pre, i)
