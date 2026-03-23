import dolfin as df
import numpy as np
from petsc4py import PETSc
from dolfin import as_backend_type

def rec_with_reg(
    u_noised,
    mesh_params={},
    mat_params={},
    Nx_coarse = int):

    Nx = mesh_params.get("Nx")
    Ny = mesh_params.get("Ny")
    degree = mesh_params.get("degree", 3)  
    mesh = df.UnitSquareMesh(Nx, Ny)
    dx = df.Measure("dx", domain=mesh)   
    tol=1e-15

    rho = mat_params.get("rho", 1)
    omega = mat_params.get("omega", 1)
    mesh_coarse = df.UnitSquareMesh(Nx_coarse, Nx_coarse)
    dx_coarse = df.Measure("dx", domain=mesh_coarse)
    degree_coarse = degree
    V1_coarse = df.VectorFunctionSpace(mesh_coarse, 'DG', degree_coarse)

    # L2 projection of u onto the coarser mesh
    u=df.TrialFunction(V1_coarse)
    v=df.TestFunction(V1_coarse)
    a_proj = df.dot(u,v)*dx_coarse 
    l_proj = df.dot(u_noised,v)*dx_coarse 
    u_proj= df.Function(V1_coarse)
    df.solve(a_proj== l_proj, u_proj,  solver_parameters={'linear_solver': 'mumps'})
    
    # Compute S_11delta, S_12delta and S_22delta
    V1_DG= df.VectorFunctionSpace(mesh,  'DG',degree)
    epsilon_udelta =  0.5*( df.nabla_grad(df.interpolate(u_proj,V1_DG)) + df.nabla_grad(df.interpolate(u_proj,V1_DG) ).T)        
    S_11delta = epsilon_udelta [0, 0]
    S_22delta = epsilon_udelta [1, 1]
    S_12delta = epsilon_udelta [0, 1]

    #Projection onto continuous FE space
    V1_G= df.FunctionSpace(mesh,  'CG',degree-1)
    S_11delta_cont = df.project(S_11delta,V1_G )
    S_22delta_cont = df.project(S_22delta, V1_G )
    S_12delta_cont = df.project(S_12delta, V1_G ) 

    #Compute the source term G
    V_S = df.FunctionSpace(mesh, 'CG',degree)  
    u1= df.interpolate(u_proj,V1_DG)[0]
    u2= df.interpolate(u_proj,V1_DG)[1]
    F=- rho * omega**2* (df.grad(u1 )[1] - df.grad( u2 )[0] )  - 2*df.grad(df.grad(S_11delta - S_22delta)[1])[0] - 2*df.grad(df.grad(S_12delta)[1])[1] + 2*df.grad(df.grad(S_12delta)[0])[0]
    def boundary(x, on_boundary):
        return on_boundary and x[0]<tol 
    bc = df.DirichletBC(V_S, df.Constant(0), boundary)
    w1 = df.TrialFunction(V_S)
    w1_test = df.TestFunction(V_S)
    a2 = df.grad(w1)[0]*w1_test*dx   
    l2 = F* w1_test*dx
    prim = df.Function(V_S)
    df.solve(a2 == l2, prim, bc)

    V_rot = df.VectorFunctionSpace(mesh, 'CG', degree-1)
    V_S1 = df.FunctionSpace(mesh, 'CG',degree-1) 
    mu1 = df.TrialFunctions(V_rot)
    nu1 = df.TestFunctions(V_rot)

    a1=df.dot(df.grad(S_12delta_cont*mu1[0])[1] + 0.5*df.grad((S_11delta_cont-S_22delta_cont)/S_12delta_cont)[0] * mu1[0]*S_12delta_cont + 0.5*(S_11delta_cont-S_22delta_cont)/S_12delta_cont*df.grad(mu1[0]*S_12delta_cont)[0] - df.grad(mu1[1]*S_12delta_cont)[0] ,nu1[0] )*dx  +  df.dot(df.grad(S_12delta_cont*mu1[1])[1] ,nu1[1] )*dx - 0.5*df.dot((S_11delta_cont-S_22delta_cont)*mu1[0],df.grad(nu1[0])[0])*dx  + df.dot(S_12delta_cont*mu1[0], df.grad(nu1[1])[0])*dx    

    l1=df.dot (prim /2,nu1[1])* dx

    def custom_boundary1(x, on_boundary):
        return on_boundary and (x[0]<tol or x[0] > 1-tol or x[1]<tol)
    def custom_boundary2(x, on_boundary):
        return on_boundary  and x[1]<tol
    bc_mu1 = df.DirichletBC(V_rot.sub(0), df.Constant(0), custom_boundary1)    
    bc_mu2 = df.DirichletBC(V_rot.sub(1), df.Constant(0), custom_boundary2)
    bc = [bc_mu1, bc_mu2]

    #Normal equation
    A = df.assemble(a1)
    b = df.assemble(l1)
    [bc.apply(A, b) for bc in bc]

    mu1 = df.Function(V_rot )
    [bc.apply(A, b) for bc in bc]

    # A^T A et A^T b
    A_petsc = as_backend_type(A).mat()
    b_petsc = as_backend_type(b).vec()
    AT_petsc = PETSc.Mat().createTranspose(A_petsc)
    ATA_petsc = AT_petsc * A_petsc
    AT_b_petsc = AT_petsc * b_petsc

    # Solution vector
    mu1_petsc = b_petsc.duplicate()

    # Configure the PETSc solver
    ksp = PETSc.KSP().create()
    ksp.setOperators(ATA_petsc)
    ksp.setFromOptions()
    ksp.setType('cg')
    ksp.getPC().setType('none')

    # Solve the linear system
    ksp.solve(AT_b_petsc, mu1_petsc)

    converged_reason = ksp.getConvergedReason()

    if converged_reason > 0:
        print("The solver has converged successfully.")
    else:
        print("Solver diverged due to high residual, but let's check the result.")
        print(f"Final residual norm: {ksp.getResidualNorm()}")
    
    mu1 = df.Function(V_rot)
    mu1.vector().set_local(mu1_petsc.getArray())

    mu = mat_params.get("mu", 2)
    mu_function = df.project(mu, V_S1)

    mu_function_rec  = df.project(mu1[0]+ 1 , V_S1) 
    error_l2 = df.sqrt(df.assemble((mu_function - mu_function_rec)**2 * dx)) / df.sqrt(df.assemble(mu_function**2 * dx))
    print('L2-norm error for the noisy approach with regularization', error_l2) 

    return mu_function_rec, error_l2

    