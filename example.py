from eqsystem import EqSystem, EqSystemGraph, PyScriptGenerator


S = EqSystem()
# Начальное приближение можно задать для всех переменных
S.add_var('x0', 0.5)
S.add_var('x1', 0.5)
S.add_var('x2', 0.5)
S.add_var('x3', 0.5)

eq0 = {'x0': (True, f"x0 = x1**2 + x3 * x2**3"),
       'x1': (True, f"x1 = (x0 - x3 * x2**3)**0.5"),
       'x2': (True, f"x2 = ((x0 - x1**2)/x3)**(1/3)"),
       'x3': (False, f"x1**2 + x * x2**3 - x0")}
eq1 = {'x0': (True, f"x0 = x2**0.5 - x1"),
       'x1': (True, f"x1 = x2**0.5 - x0"),
       'x2': (True, f"x2 = (x0 + x1)**2")}
eq2 = {'x0': (True, f"x0 = x2"),
       'x2': (True, f"x2 = x0")}
S.add_std_eq('eq0', eq0)
S.add_std_eq('eq1', eq1)
S.add_std_eq('eq2', eq2)

S.set_desired_x({'x3'})
S.set_desired_y({'x1'})

G = EqSystemGraph(S)
G.solve_matching_maxflow_ortools()
G.scc()
G.mfas()
gen = PyScriptGenerator("out.py", "solve")
G.generate_code(gen)