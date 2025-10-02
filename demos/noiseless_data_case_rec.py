import dolfin as df
import numpy as np

def noiseless_data_rec(
        u, 
        mesh_params={},
        mat_params={},
        load_params={}):

    Nx = mesh_params.get("Nx")
    Ny = mesh_params.get("Ny")
    degree = mesh_params.get("degree", 3)
    mesh = df.UnitSquareMesh(Nx, Ny)
    dx = df.Measure("dx", domain=mesh)
    V_S2 = df.FunctionSpace(mesh, 'DG',degree-1)  
    def epsilon(u):
        return 0.5*(df.nabla_grad(u) + df.nabla_grad(u).T) 
    V = df.VectorElement('CG', mesh.ufl_cell(), degree=degree)
    V1 = df.FunctionSpace(mesh, V)
    u = df.interpolate(u, V1)
    S_11 = epsilon(u)[0, 0]
    S_22 = epsilon(u)[1, 1]
    S_12 = epsilon(u)[0, 1]
    S = (S_11-S_22)/S_12

    mu = mat_params.get("mu", 2)
    rho = mat_params.get("rho", 1)
    omega = mat_params.get("omega", 1)
    f = load_params.get("f")

    f1=f[0]
    f2=f[1]
    u1=u[0]
    u2=u[1]
    V_S = df.FunctionSpace(mesh, 'CG',degree)  

    ### Compute the antiderivative of $g$ with respect to $x_2$, denoted by $G$ in the article.
    F = -rho * omega**2 * df.grad(u1 + f1)[1] + rho * omega**2 * df.grad(u2 + f2)[0] - 2*df.grad(df.grad(S_11 - S_22)[1])[0] - 2*df.grad(df.grad(S_12)[1])[1] + 2*df.grad(df.grad(S_12)[0])[0]
    
    tol = 1e-15
    def boundary(x, on_boundary):
         return on_boundary and x[1]<tol 
    bc = df.DirichletBC(V_S, df.Constant(0), boundary)
    w1 = df.TrialFunction(V_S)
    w1_test = df.TestFunction(V_S) 
    a2 = df.grad(w1)[1] * w1_test * dx   
    l2 = F * w1_test * dx
    prim = df.Function(V_S)
    df.solve(a2 == l2, prim, bc, solver_parameters={'linear_solver':'mumps'})
    
    V_rot = df.VectorFunctionSpace(mesh, 'CG', degree-1)
    V_S1 = df.FunctionSpace(mesh, 'CG', degree-1) 
    mu1 = df.TrialFunctions(V_rot)
    nu1 = df.TestFunctions(V_rot)

    # Normal equation
    a1 = df.dot(df.grad(mu1[0]*(S_11-S_22))[0] + df.grad(mu1[0]*S_12)[1] - df.grad(mu1[1]*S_12)[0], df.grad(nu1[0]*(S_11-S_22))[0] + df.grad(nu1[0]*S_12)[1] - df.grad(nu1[1]*S_12)[0]) * dx  + df.dot(df.grad(mu1[1]*S_12)[1] - df.grad(mu1[0]*S_12)[0], df.grad(nu1[1]*S_12)[1] - df.grad(nu1[0]*S_12)[0]) * dx
    l1 = df.dot(prim/2, df.grad(nu1[0]*(S_11-S_22))[0] + df.grad(nu1[0]*S_12)[1] - df.grad(nu1[1]*S_12)[0]) * dx

    def custom_boundary1(x, on_boundary):
        return on_boundary and (x[0]<tol or x[0]>1-tol or x[1]<tol)
    def custom_boundary2(x, on_boundary):
        return on_boundary  and x[1]<tol 
    bc_mu1 = df.DirichletBC(V_rot.sub(0), df.Constant(0), custom_boundary1) 
    bc_mu2 = df.DirichletBC(V_rot.sub(1), df.Constant(0), custom_boundary2)
    bc = [bc_mu1, bc_mu2]

    mu1 = df.Function(V_rot)
    df.solve(a1 == l1, mu1, bc, solver_parameters={'linear_solver':'mumps'})

    mu_function_rec = df.project(mu1[0]+1, V_S1)  

    mu_function = df.project(mu, V_S1)
    norm1_l2 = df.sqrt(df.assemble((mu_function - mu_function_rec)**2 * dx)) / df.sqrt(df.assemble(mu_function**2 * dx))
    print('For noiseless data, the L2-norm for the reconstruction:', norm1_l2)

    return mu_function_rec, norm1_l2, prim
