/*
	File:			Solve.h

	Function:		Contains routines for solving a system of linear equations.
					Includes the overrelaxation (a more general version of 
					Gauss Seidel) and conjugate gradient methods, for both
					normal and sparse matrices.

	Author(s):		Andrew Willmott

	Copyright:		(c) 1995-2000, Andrew Willmott
 */

#ifndef __Solve__
#define __Solve__

#include "vl/VL.h"
#include "vl/Mat.h"
#include "vl/SparseMat.h"

TMReal SolveOverRelax(const TMat &A, TVec &x, const TVec &b,
			TMReal epsilon, TMReal omega = 1.0, Int *steps = 0);
TMReal SolveOverRelax(const TSparseMat &A, TVec &x, const TVec &b, 
			TMReal epsilon, TMReal omega = 1.0, Int *steps = 0);

TMReal SolveConjGrad(const TMat &A, TVec &x, const TVec &b,
			TMReal epsilon, Int *steps = 0); 
TMReal SolveConjGrad(const TSparseMat &A, TVec &x, const TVec &b,
			TMReal epsilon, Int *steps = 0);

TMReal SolveConjGrad(
				const TSparseMat	&A,
				TVec				&x,
				const TVec			&b,
				TMReal				epsilon,
				Int					*steps
			)
{
	Int		i, iMax;
	TMReal	alpha, beta, rSqrLen, rSqrLenOld, u;
	TVec	r(b.Elts());		// Residual vector, b - Ax 
	TVec	d(b.Elts());		// Direction to move x at each step 
	TVec	t(b.Elts());

	Assert(A.IsSquare(), "(SolveConjGrad) Matrix not square");
	r = b;
	r -= A * x;				
	rSqrLen = sqrlen(r);
	d = r;						
	i = 0;
	if (steps)
		iMax = *steps;		
	else
		iMax = VL_MAX_STEPS;
		
	if (rSqrLen > epsilon)				// If we haven't already converged...
		while (i < iMax)
		{	
			i++;
			t = A * d;		
			u = dot(d, t);
			
			if (u == 0.0)
			{
				_Warning("(SolveConjGrad) d'Ad = 0");
				break;
			}
			
			alpha = rSqrLen / u;		// How far should we go?
			x += alpha * d;				// Take a step along direction d
			r -= alpha * t; 
			rSqrLenOld = rSqrLen;
			rSqrLen = sqrlen(r); 
			
			if (rSqrLen <= epsilon)
				break;					// Converged! Let's get out of here
			
			beta = rSqrLen / rSqrLenOld;
			d = r + beta * d;			//	Change direction
		}

	if (steps)
		*steps = i;
	
	return(rSqrLen);
}

#endif
