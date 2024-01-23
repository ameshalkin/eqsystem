import igraph as ig
import numpy as np
from ortools.graph.python import max_flow
import time


class EqSystem:
    def __init__(self):
        self.var_ids = {}
        self.vars_used = []
        self.X0 = {}

        self.eq_ids = {}
        self.eqs = []

        self.desired_x = set()
        self.desired_y = set()

    def add_var(self, var_name, x0 = None):
        if var_name in self.var_ids.keys():
            raise ValueError(f"Variable name '{var_name}' already exists")
        if not var_name.isidentifier():
            raise ValueError(f"Variable name '{var_name}' not valid")
        self.var_ids[var_name] = len(self.var_ids)
        self.vars_used.append(False)
        if x0:
            self.X0[var_name] = x0

    def add_std_eq(self, eq_name, equation):
        if eq_name in self.eq_ids:
            raise ValueError(f"Equation name '{eq_name}' already exists")
        if not eq_name.isidentifier():
            raise ValueError(f"Variable name '{eq_name}' not valid")
        for var in equation.keys():
            if var not in self.var_ids.keys():
                raise ValueError(f"Unknown variable name: '{var}'")
            if not self.vars_used[self.var_ids[var]]:
                self.vars_used[self.var_ids[var]] = True
        self.eq_ids[eq_name] = len(self.eq_ids)
        self.eqs.append(equation)

    def set_desired_x(self, desired_x: set):
        for var in desired_x:
            if var not in self.var_ids:
                raise ValueError(f"Unknown variable name: '{var}'")
        self.desired_x = desired_x

    def set_desired_y(self, desired_y: set):
        for var in desired_y:
            if var not in self.var_ids:
                raise ValueError(f"Unknown variable name: '{var}'")
        self.desired_y = desired_y
    
    def get_expr(self, eq_name):
        return self.eqs[self.eq_ids[eq_name]]

    
class EqSystemGraph:
    def __init__(self, eqsystem):
        self.eqsystem = eqsystem
        self.m = len(eqsystem.eq_ids)
        self.n = len(eqsystem.var_ids)
        num_vars_used = sum(eqsystem.vars_used)
        if num_vars_used != len(eqsystem.var_ids): 
            raise ValueError(f"Variables used: {num_vars_used}, "
                             f"number of variables: {len(eqsystem.var_ids)}")
        edge_list = []
        for _, id in eqsystem.eq_ids.items():
            edge_list += [(id, self.m + eqsystem.var_ids[var]) for var in eqsystem.eqs[id].keys()]
        self.g = ig.Graph(self.m+self.n, edge_list, directed=True)
        self.g.vs["name"] = list(eqsystem.eq_ids.keys()) + list(eqsystem.var_ids.keys())
        self.g.vs["type"] = ["eq"] * self.m + ["var"] * self.n
        self.desired_x = set(map(lambda x: self.m + eqsystem.var_ids[x], eqsystem.desired_x))
        self.desired_y = set(map(lambda y: self.m + eqsystem.var_ids[y], eqsystem.desired_y))
        self.eqs_range = range(self.m)
        self.vars_range = range(self.m, self.m + self.n)

    def solve_matching_maxflow_ortools(self):

        #ts1 = time.time()

        smf = max_flow.SimpleMaxFlow()
        v = self.g.vcount()
        start_nodes = []
        end_nodes = []
        for u, w in self.g.get_edgelist():
            start_nodes.append(u)
            end_nodes.append(w)

        start_nodes += [v for i in self.eqs_range]
        end_nodes += [i for i in self.eqs_range]

        for var in self.vars_range:
            if var in self.desired_y:
                start_nodes.append(var)
                end_nodes.append(v + 2)
            elif var not in self.desired_x:
                start_nodes.append(var)
                end_nodes.append(v + 1)
        
        start_nodes.append(v+1)
        end_nodes.append(v+2)

        capacities = np.ones(len(start_nodes))
        capacities[-1] = self.m - len(self.desired_y)
        start_nodes = np.array(start_nodes)
        end_nodes = np.array(end_nodes)
        all_arcs = smf.add_arcs_with_capacity(start_nodes, end_nodes, capacities)

        status = smf.solve(v, v+2)
        
        if status != smf.OPTIMAL:
            print("There was an issue with the max flow input.")
            print(f"Status: {status}")
            exit(1)
        #print("Max flow:", smf.optimal_flow())
        if smf.optimal_flow() != self.m:
            raise ValueError(f"Equation-variable matching of size {self.m} does not exist\n" \
                             f"Max flow: {smf.optimal_flow()}") 

        solution_flows = smf.flows(all_arcs)

        zero_edges = np.where(solution_flows[:self.g.ecount()] == 0)[0]
        self.g.reverse_edges(zero_edges.tolist())

        edges_to_add = []
        vs_to_remove = []
        for eq in self.eqs_range:
            next_var = self.g.successors(eq)[0]
            next_var_name = self.g.vs[next_var]["name"]
            self.g.vs[eq]["name"] += f"-{next_var_name}"
            edges_to_add += [(eq, next_eq) for next_eq in self.g.successors(next_var)]
            vs_to_remove.append(next_var)
        self.g.add_edges(edges_to_add)
        self.g.delete_vertices(vs_to_remove)

        #ts2 = time.time()
        #print(f"Maxflow OR-Tools time: {ts2 - ts1}")

    def scc(self):
        #ts1 = time.time()
        self.sccs = self.g.connected_components(mode='strong')
        #ts2 = time.time()
        #print(f"SCC time: {ts2 - ts1}")

    def mfas(self):
        #ts1 = time.time()

        self.no_cycle_components = {}
        for i, cluster in enumerate(self.sccs):
            if len(cluster) > 1:
                component = self.sccs.subgraph(i)

                mfas = component.feedback_arc_set(method="eades")
                vars_to_init = []
                for edge in mfas:
                    id = component.get_edgelist()[edge][0]
                    vars_to_init.append(component.vs[id]["name"].split('-')[1])

                component.delete_edges(mfas)
                ordering = component.topological_sorting(mode='out')
                eqs_order = [component.vs[i]["name"] for i in ordering]

                self.no_cycle_components[i] = (vars_to_init, eqs_order)

        #ts2 = time.time()
        #print(f"MFAS: {ts2 - ts1}")
    
    def generate_code(self, gen):
        
        input_idx = 0
        self.input_ids = {}
        self.input_names = []
        self.output_ids = {}
        self.output_names = []

        for i, cluster in enumerate(self.sccs):
            if len(cluster) > 1:
                vars_to_init, eqs_order = self.no_cycle_components[i]

                init_vars = {}
                for var_name in vars_to_init:
                    x0 = self.eqsystem.X0[var_name]
                    if x0 is None: raise ValueError(f"No initial value for '{var_name}'")
                    init_vars[var_name] = x0

                exprs = []
                for eq in eqs_order:
                    eq_name, var_name = eq.split('-')
                    inverse_flag, expr = self.eqsystem.get_expr(eq_name)[var_name]
                    if not inverse_flag:
                        x0 = self.eqsystem.X0[var_name]
                        if x0 is None: raise ValueError(f"No initial value for '{var_name}'")
                    else:
                        x0 = None
                    exprs.append((var_name, x0, eq_name, inverse_flag, expr))

                gen.add_cycle(init_vars, exprs, indent=0)

            elif self.g.vs[cluster[0]]["type"] == 'var': 
                var_name = self.g.vs[cluster[0]]["name"]
                gen.add_variable_assign(var_name, input_idx, indent=0)
                self.input_ids[var_name] = input_idx
                self.input_names.append(var_name)
                input_idx += 1

            elif self.g.vs[cluster[0]]["type"] == 'eq':
                eq_name, var_name = self.g.vs[cluster[0]]["name"].split('-')
                inverse_flag, expr = self.eqsystem.get_expr(eq_name)[var_name]
                if inverse_flag:
                    gen.add_simple_equation(expr, indent=0)
                else:
                    x0 = self.eqsystem.X0[var_name]
                    if x0 is None: raise ValueError(f"No initial value for '{var_name}'")
                    gen.add_hard_equation(eq_name, var_name, x0, 
                                             expr, indent=0)
                    
        for i, eq in enumerate(self.g.vs[:self.m]):
            _, var_name = eq["name"].split('-')
            self.output_names.append(var_name)
            self.output_ids[var_name] = i
        gen.add_output_assign(self.output_names, indent=0)

        gen.write_to_file()

class PyScriptGenerator:
    def __init__(self, filename, name):
        self.filename = filename
        self.head = 'import numpy as np\nfrom scipy.optimize import fsolve\n'
        if not name.isidentifier(): raise ValueError(f"Function name '{name}' is not valid")
        self.name = name
        self.func_def = f'def {name}(x):\n'
        self.defs = ''
        self.func_text = ''
        self.func_indent = 1
        self.rel_eps = 0.0001
        self.abs_eps = 1e-8
    
    def get_indent(self, n = 0):
        return (self.func_indent + n) * '    '
    
    def add_variable_assign(self, var_name, idx, indent=0):
        self.func_text += self.get_indent(indent) + f"{var_name} = x[{idx}]" + '\n'

    def add_variable_init(self, var_name, x0, indent=0):
        self.func_text += self.get_indent(indent) + f"{var_name} = {x0}" + '\n'

    def add_cycle(self, init_vars, exprs, indent=0):
        for var_name, x0 in init_vars.items():
            self.add_variable_init(var_name, x0, indent=indent)

        self.func_text += self.get_indent() + f"while True:" + '\n'

        for var_name, x0 in init_vars.items():
            self.add_variable_init("temp_"+var_name, var_name, indent=indent+1)

        for var_name, x0, eq_name, inverse_flag, expr in exprs:
            if inverse_flag:
                self.add_simple_equation(expr, indent=indent+1)
            else:
                self.add_hard_equation(eq_name, var_name, x0, expr, indent=indent+1)

        for var_name in init_vars.keys():
            self.func_text += self.get_indent(indent+1) + \
                              f"if abs(({var_name} - temp_{var_name})" \
                              f" / (temp_{var_name} + {self.abs_eps})) > {self.rel_eps}: continue\n"
        self.func_text += self.get_indent(1) + "break\n"


    def add_simple_equation(self, formula, indent=0):
        self.func_text += self.get_indent(indent) + formula + '\n'
    
    def add_hard_equation(self, eq_name, var_name, x0, formula, indent=0):
        if eq_name == self.name: raise ValueError(f"Function names coincide: '{self.name}'")
        self.defs += self.get_indent(0) + f"def {eq_name}(x):\n" +\
                     self.get_indent(0) + f"    return {formula}\n"

        self.func_text += self.get_indent(indent) + \
                          f"{var_name} = fsolve({eq_name}, {x0}, xtol = {self.rel_eps})[0]\n"
    
    def add_output_assign(self, out_names, indent=0):
        self.func_text += self.get_indent(indent) + f"y = [0.0] * {len(out_names)}" + '\n'
        for i, var_name in enumerate(out_names):
            self.func_text += self.get_indent(indent) + f"y[{i}] = {var_name}" + '\n'
        self.func_text += self.get_indent(indent) + "return y"

    def write_to_file(self):
        with open(self.filename, "w") as file:
            file.write(self.head + '\n')
            file.write(self.func_def)
            file.write(self.defs + '\n')
            file.write(self.func_text)
   