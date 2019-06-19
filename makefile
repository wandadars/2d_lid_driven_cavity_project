all:
	gfortran -fbounds-check -frepack-arrays mod_parameters.f95 cg_method.f95 CalcConvergence.f95 SetInitCond.f95 TimeStep.f95 UpdateStars.f95 WriteSolution.f95 StarSolver.f95 ReadInput.f95 PressPoisson.f95 UpdateVelocity.f95 main.f95 -o 2DNavierStokes
