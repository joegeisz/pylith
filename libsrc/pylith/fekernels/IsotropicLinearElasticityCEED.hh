
#ifndef isotropiclinearelasticityceed_h
#define isotropiclinearelasticityceed_h



/*
============================================================================================================================================
structs used to pass in constant parameters
============================================================================================================================================
*/

struct CoordinatesContext {

	CeedInt dim;
	///...

}///end CoordinatesContext definition

struct ProblemContext {

	//dimension
	CeedInt dim;

	//Bulk Modulus
	CeedScalar BM;

	//Shear Modulus
	CeedScalar SM;
	///...

}///end ProblemContext definition



/*
============================================================================================================================================
functions for making Qfunctions
============================================================================================================================================
*/


/*
/////////////////////////////////////////////////////////////////////////////////////////

Setup calculates the inverse of the jacobian in order to translate between local and 
global coordinates. Works for 2d or 3d cases. 

Input:
	J = jacobian
	w = weights

Output:
	qdata = weight * determinant of jacobian * inverse of jacobian

/////////////////////////////////////////////////////////////////////////////////////////
*/
static int Setup(void *ctx, CeedInt Q, const CeedScalar *const *in, CeedScalar *const *out) {

	//name the parameters from ctx
	struct CoordinatesContext *context = (struct CoordinatesContext *)ctx;
	//...

	//name the input vectors
	const CeedScalar *J = in[0];
	const CeedScalar *w = in[1];
	//...

	//name the output vectors
	CeedScalar *qdata = out[0];
	//...

	//loop over quadrature points, depending on problem dimensionality
	if(context.dim == 2){
		//------------------------------------
		//Calculate 2D inverse
		//------------------------------------
		for(int i=0; i<Q; i++){
			
			//Setup
			const CeedScalar a = J[i+0*Q];
			const CeedScalar b = J[i+1*Q];
			const CeedScalar c = J[i+2*Q];
			const CeedScalar d = J[i+3*Q];
			//const CeedScalar detJ = J[i+Q*0]*J[i+Q*3]-J[i+Q*2]*J[i+Q*1];

			// Qdata
			// This is the inverse of the Jacobian times its determinant times the weight:
    		qdata[i+ 0*Q] = w[i] * d;
			qdata[i+ 1*Q] = w[i] * -b;
			qdata[i+ 2*Q] = w[i] * -c;
			qdata[i+ 3*Q] = w[i] * a;
    	}//end 2d for loop
	}//end 2d case

	else if(context.dim == 3){
		//------------------------------------
		// Calculate 3d inverse
		//------------------------------------
		for(int i=0; i<Q; i++){

		    // Setup
		    const CeedScalar J11 = J[i+Q*0];
		    const CeedScalar J21 = J[i+Q*1];
		    const CeedScalar J31 = J[i+Q*2];
		    const CeedScalar J12 = J[i+Q*3];
		    const CeedScalar J22 = J[i+Q*4];
		    const CeedScalar J32 = J[i+Q*5];
		    const CeedScalar J13 = J[i+Q*6];
		    const CeedScalar J23 = J[i+Q*7];
		    const CeedScalar J33 = J[i+Q*8];
		    const CeedScalar A11 = J22*J33 - J23*J32;
		    const CeedScalar A12 = J13*J32 - J12*J33;
		    const CeedScalar A13 = J12*J23 - J13*J22;
		    const CeedScalar A21 = J23*J31 - J21*J33;
		    const CeedScalar A22 = J11*J33 - J13*J31;
		    const CeedScalar A23 = J13*J21 - J11*J23;
		    const CeedScalar A31 = J21*J32 - J22*J31;
		    const CeedScalar A32 = J12*J31 - J11*J32;
		    const CeedScalar A33 = J11*J22 - J12*J21;
		    //const CeedScalar detJ = J11*A11 + J21*A12 + J31*A13;


    		// Qdata
			// This is the inverse of the Jacobian times its determinant times the weight:
			qdata[i+ 0*Q] = w[i] * A11;
			qdata[i+ 1*Q] = w[i] * A12;
			qdata[i+ 2*Q] = w[i] * A13;
			qdata[i+ 3*Q] = w[i] * A21;
			qdata[i+ 4*Q] = w[i] * A22;
			qdata[i+ 5*Q] = w[i] * A23;
			qdata[i+ 6*Q] = w[i] * A31;
			qdata[i+ 7*Q] = w[i] * A32;
			qdata[i+ 8*Q] = w[i] * A33;
		}//end 3d for loop
	}//end 3d case

	else{
		printf("Dimension not valid\n")
	}
	//end loop over quadrature points

	return 0;

}///end Setup definition



/*
/////////////////////////////////////////////////////////////////////////////////////////
ApplyMass applies the mass matrix to the vector v. Works for 2D or 3D. 

Inputs: 
	u = displacement vector
	w = weights

Output:
	v = w.v

/////////////////////////////////////////////////////////////////////////////////////////
*/
static int ApplyMass(void *ctx, CeedInt Q, const CeedScalar *const *in, CeedScalar *const *out) {

	//name the parameters from ctx
	struct CoordinatesContext *context = (struct CoordinatesContext *)ctx;
	//...

	//name the input vectors
	const CeedScalar *u = in[0];
	const CeedScalar *w = in[1];
	//...

	//name the output vectors
	CeedScalar *v = out[0];
	//...

	//loop over quadrature points
	if(context.dim == 2){
		for(int i=0; i<Q; i++){
			/*
			//multiply weights
			*/

	    	v[i+0*Q] = w[i+0*Q] * u[i+0*Q];
	    	v[i+1*Q] = w[i+1*Q] * u[i+1*Q];

		}//end loop over quadrature points
	}
	else if(context.dim = 3){
		for(int i=0; i<Q; i++){
			/*
			//multiply weights
			*/

	    	v[i+0*Q] = w[i+0*Q] * u[i+0*Q];
	    	v[i+1*Q] = w[i+1*Q] * u[i+1*Q];
	    	v[i+2*Q] = w[i+2*Q] * u[i+2*Q];

		}//end loop over quadrature points

	}

	return 0;

}///end ApplyMass definition


/////////////////////////////////////////////////////////////////////////////////////////
/*
IsotropicLinearElasticity implements the physics of the Linear Elasticity equations either in 
full 3d or in the 2d case of plane strain. 

Governing equation:

	s_ij,j + f_i = rho (d2u/dt2) in volume

Weak Form:

		I[ -s_ij*vi,j ]dV 
	+ 	I[ f_i*v_i ]dV 
	+	I[ -rho*(d2u_i/dt2)*v_i ]dV 
	+ 	I[ T_i*v_i ]dS 
	= 0

	2d term (body force/gravity) is delt with in other parts of code,
	3rd term (time dependent) is delt with by time stepping,
	4th term (surface) is delt with by a different integrator, 
	so we only worry about the first term.



Input vector 2D plane strain:

	displacement u = (ux, uy)

	displacment gradient du = (	dux/dx, duy/dx, 
								dux/dy, duy/dy 	)

	coordinate transform matrix wBJ = inverse of Jacobian from setup


Output vector 2D plane strain:

	no terms multiplied by test function only:
		v = (0, 0) 

	terms multiplied by gradient of test function (but adjusted by wBJ to account for test coordinates):
		dv = 	(s11 s12)
				(s21 s22)

Input vector 3D:

	displacement u = (ux, uy, uz)

	displacement gradient du = (dux/dx, duy/dx, duz/dx
								dux/dy, duy/dy, duz/dy
								dux/dz, duy/dz, duz/dz)

	coordinate transform matrix wBJ = inverse of Jacobian from setup 

Output vector 3D:

	no terms multiplied by test function only:
		v = (0, 0, 0) 

	terms multiplied by gradient of test function (but adjusted by wBJ to account for test coordinates):
		dv = 	(s11 s12 s13)
				(s21 s22 s23)
				(s21 s22 s23)

Math:

	Symmetrize the du matrix for the strain tensor:

		Strain Tensor E = 1/2 (du + du^T)

	Then find stress tensor from strain tensor:

		sigma_ij = C_ijkl E_kl 

		in notation from Classical Mechanics (Taylor):

			sigma = ae1 + bEp

			where:
				a = 3*K
				b = 2*G
				K = Bulk Modulus
				G = Shear Modulus
				e = (1/3)*tr(E) = "Spherical Part of strain"
				Ep = E - e1 = "Deviatoric Part of strain"

		--OR--

		in voight notation, using symmetries:

			sigma_i = D_ij * E_j

			where:
				sigma 	= (s11, s22, s33, s23, s13, s12)T
				E 		= (E11, E22, E33, E23, E13, E12)T
				D_ij = 	[	2*mu+lambda		lambda		lambda		0	0	0	0	]
						[	lambda			2*mu+lambda	lambda		0	0	0	0	]
						[	lambda			lambda		2*mu+lambda	0	0	0	0	]
						[	0				0			0			mu 	0	0	0	]
						[	0				0			0			0 	mu	0	0	]
						[	0				0			0			0 	0	mu	0	]
						[	0				0			0			0 	0	0	mu 	]
				mu = G
				lambda = K - (2/3)*G






*/
/////////////////////////////////////////////////////////////////////////////////////////

static int IsotropicLinearElasticity(void *ctx, CeedInt Q, const CeedScalar *const *in, CeedScalar *const *out) {

	//name the parameters from ctx
	struct ProblemContext *context = (struct ProblemContext *)ctx;
	//...

	//name the input vectors
	const CeedScalar *disp = in[0];
	const CeedScalar *graddisp = in[1];
	const CeedScalar *qdata = in[2];
	//...

	//name the output vectors
	CeedScalar *v = out[0];
	CeedScalar *dv = out[0];
	
	//...

	//loop over quadrature points


	//------------------------------------
	//2D Linear Elastic equations, plane strain
	//------------------------------------
	if(context.dim==2){
			for(int i = 0; i<Q, i++){
			//Linear Elastic equations for 2d Plane Strain problem
			const CeedScalar u[2] = { disp[i+0*Q], disp[i+1*Q] };

			const CeedScalar du[4] = { 	graddisp[i+0*Q], 
										graddisp[i+1*Q], 
										graddisp[i+2*Q], 
										graddisp[i+3*Q] };

			//strain tensor
			const CeedScalar E[3] = { 	du[0],
										du[3],
										.5*(du[1]+du[2])	};

			//Relate Stress and strain
			const CeedScalar mu = context.SM;
			const CeedScalar lam = context.BM - (2/3)*context.SM;


			const CeedScalar sigma[4] = {	(lam+2*mu)E[0] + lam*E[1],		//Sigma11
											lam*E[0] + (lam+2*mu)E[1],		//Sima22
											2*mu*E[3],						//Sigma12
											lam*(E[0]+E[1])					//Sigma33 ???
												};
			//set outputs

			//body force dealt with elsewhere, only need sigma term
			v[i+0*Q] = 0;
			v[i+1*Q] = 0;
			
			//-sigmaij*vi,j
			dv[i+0*Q] = -(sigma[0]*wBJ[0]+sigma[2]*wBJ[1]); //v11
			dv[i+1*Q] = -(sigma[2]*wBJ[0]+sigma[3]*wBJ[1]); //v12

			dv[i+3*Q] = -(sigma[0]*wBJ[2]+sigma[2]*wBJ[3]); //v21
			dv[i+4*Q] = -(sigma[2]*wBJ[2]+sigma[3]*wBJ[3]); //v22

			//??? dv[i+5*Q] = -sigma[3]; //v33


			}//end 2d foor loop
	}//end 2d case

	//------------------------------------
	//Full 3d Linear Elastic equations
	//------------------------------------
	else if(context.dim == 3){
		for(int i=0; i<Q; i++){
		
			//Displacement
			const CeedScalar u[3] = { disp[i+0*Q], disp[i+1*Q], dis[i+2*Q] };

			//Grad Displacement
			const CeedScalar du[9] = { 	graddisp[i+0*Q], graddisp[i+1*Q], graddisp[i+2*Q],
										graddisp[i+3*Q], graddisp[i+4*Q], graddisp[i+5*Q],
										graddisp[i+6*Q], graddisp[i+7*Q], graddisp[i+8*Q] };

			//jacobian/coordinate transform data
			const CeedScalar wBJ[9] = {	qdata[i+0*Q], qdata[i+1*Q], qdata[i+2*Q],
										qdata[i+3*Q], qdata[i+4*Q], qdata[i+5*Q],
										qdata[i+6*Q], qdata[i+7*Q], qdata[i+8*Q] };

			//------------------------------------
			//Full Symmetric 3x3 Matrix Notation (less efficient, more readable, follows notation from Taylor's Classical Mechanics)		
			//------------------------------------
			if(1){

				
				//Strain Tensor 
				const CeedScalar E[9] = {		du[0], 			.5*(du[1]+du[3]), 	.5*(du[2]+du[6]),
											.5*(du[1]+du[3]), 		du[4], 			.5*(du[5]+du[7]),
											.5*(du[2]+du[6]), 	.5*(du[5]+du[7]), 		du[8] 			};

				//Spherical Part
				const CeedScalar e = (1/3)*(E[0]+E[4]+E[8]);

				//Deviatoric Part (not necessary for 2nd version of stress tensor)
				/*
				const CeedScalar Ep[9] = { 	E[0]-e, E[1], E[2],
											E[3], E[4]-e, E[5],
											E[6], E[7], E[8]-e,	};
				*/


				//Relate Stress and Strain Tensors

				const CeedScalar a = 3*context.BM;
				const CeedScalar b = 2*context.SM;


				//Stress Tensor
				/*
				const CeedScalar Sigma[9] = {	a*e + b*Ep[0], 		b*Ep[1], 				b*Ep[2],
												b*Ep[3], 			a*e + b*Ep[4], 			b*Ep[5],
												b*Ep[6], 			b*Ep[7],				a*e + b*Ep[8] };
				*/

				const CeedScalar Sigma[9] = {	(a-b)*e + b*E[0], 	b*E[1], 			b*E[2],
												b*E[3], 			(a-b)*e + b*E[4], 	b*E[5],
												b*E[6], 			b*E[7],				(a-b)*e + b*E[8] };

				//set outputs

				//body force dealt with elsewhere, only need sigma term
				v[i+0*Q] = 0;
				v[i+1*Q] = 0;
				v[i+2*Q] = 0;
			
				//-sigmaij*vi,j
				dv[i+0*Q] = -(sigma[0]*wBJ[0]+sigma[3]*wBJ[1]+sigma[6]*wBJ[2]);
				dv[i+1*Q] = -(sigma[1]*wBJ[0]+sigma[4]*wBJ[1]+sigma[7]*wBJ[2]);
				dv[i+2*Q] = -(sigma[2]*wBJ[0]+sigma[5]*wBJ[1]+sigma[8]*wBJ[2]);

				dv[i+3*Q] = -(sigma[0]*wBJ[3]+sigma[3]*wBJ[4]+sigma[6]*wBJ[5]);
				dv[i+4*Q] = -(sigma[1]*wBJ[3]+sigma[4]*wBJ[4]+sigma[7]*wBJ[5]);
				dv[i+5*Q] = -(sigma[2]*wBJ[3]+sigma[5]*wBJ[4]+sigma[8]*wBJ[5]);

				dv[i+6*Q] = -(sigma[0]*wBJ[6]+sigma[3]*wBJ[7]+sigma[6]*wBJ[8]);
				dv[i+7*Q] = -(sigma[1]*wBJ[6]+sigma[4]*wBJ[7]+sigma[7]*wBJ[8]);
				dv[i+8*Q] = -(sigma[2]*wBJ[6]+sigma[5]*wBJ[7]+sigma[8]*wBJ[8]);

			}//end full symetric 3x3 matrix notation


			//---------------------
			//Voigt Notation (utilizes symmetry to save space/calculations)
			//---------------------
			if(0){
				
				//Strain Tensor 
				const CeedScalar E[6] = {	
											du[0],
											du[4],
											du[8],
											.5*(du[5]+du[7]),
											.5*(du[2]+du[6]),
											.5*(du[1]+du[3]),
										};
								


				//Relate Stress and strain
				const CeedScalar mu = context.SM;
				const CeedScalar lam = context.BM - (2/3)*context.SM;

				//Stress Tensor
				const CeedScalar sigma[6] = { 	
												(2*mu+lam)*E[0] + lam*E[1] + lam*E[2],
												lam*E[0] + (2*mu+lam)*E[1] + lam*E[2],
												lam*E[0] + lam*E[1] + (2*mu+lam)*E[2],
												mu*E[3],
												mu*E[4],
												mu*E[5]	
											};

				//set outputs

				//body force dealt with elsewhere, only need sigma term
				v[i+0*Q] = 0;
				v[i+1*Q] = 0;
				v[i+2*Q] = 0;
			
				//-sigmaij*vi,j
				dv[i+0*Q] = -(sigma[0]*wBJ[0]+sigma[5]*wBJ[1]+sigma[4]*wBJ[2]);
				dv[i+1*Q] = -(sigma[5]*wBJ[0]+sigma[1]*wBJ[1]+sigma[3]*wBJ[2]);
				dv[i+2*Q] = -(sigma[4]*wBJ[0]+sigma[3]*wBJ[1]+sigma[2]*wBJ[2]);

				dv[i+3*Q] = -(sigma[0]*wBJ[3]+sigma[5]*wBJ[4]+sigma[4]*wBJ[5]);
				dv[i+4*Q] = -(sigma[5]*wBJ[3]+sigma[1]*wBJ[4]+sigma[3]*wBJ[5]);
				dv[i+5*Q] = -(sigma[4]*wBJ[3]+sigma[3]*wBJ[4]+sigma[2]*wBJ[5]);

				dv[i+6*Q] = -(sigma[0]*wBJ[6]+sigma[5]*wBJ[7]+sigma[4]*wBJ[8]);
				dv[i+7*Q] = -(sigma[5]*wBJ[6]+sigma[1]*wBJ[7]+sigma[3]*wBJ[8]);
				dv[i+8*Q] = -(sigma[4]*wBJ[6]+sigma[3]*wBJ[7]+sigma[2]*wBJ[8]);

			}//end voigt notation

		}//end 3d foor loop
	}//end 3d case

	else{
		printf("Dimension not valid \n");
		}//end error case


return 0;
}///end IsotropicLinearElasticity definition




#endif


