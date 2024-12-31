#pragma once
/*
* Quick and Dirty implementation of a cubic Bezier curve 
* based on https://www.codeproject.com/Articles/31859/Draw-a-Smooth-Curve-through-a-Set-of-2D-Points-wit
* 
*/
class BezierSpline
{
public:
	static void GetCurveControlPoints(double (*knots)[2], int numKnots, double (*firstControlPoints)[2], double (*secondControlPoints)[2])
	{
		if (knots == nullptr)
			return;

		int n = numKnots - 1;
		if (n < 1)
			return;

		if (n == 1)
		{ // Special case: Bezier curve should be a straight line.

			// 3P1 = 2P0 + P3
			firstControlPoints[0][0] = (2 * knots[0][0] + knots[1][0]) / 3;
			firstControlPoints[0][1] = (2 * knots[0][1] + knots[1][1]) / 3;

			// P2 = 2P1 – P0
			secondControlPoints[0][0] = 2 *
				firstControlPoints[0][0] - knots[0][0];
			secondControlPoints[0][1] = 2 *
				firstControlPoints[0][1] - knots[0][1];
			return;
		}

		// Calculate first Bezier control points
		// Right hand side vector
		double* rhs = new double [n];

		// Set right hand side X values
		for (int i = 1; i < n - 1; ++i)
			rhs[i] = 4 * knots[i][0] + 2 * knots[i + 1][0];
		rhs[0] = knots[0][0] + 2 * knots[1][0];
		rhs[n - 1] = (8 * knots[n - 1][0] + knots[n][0]) / 2.0;
		// Get first control points X-values
		double* x = GetFirstControlPoints(rhs, n);

		// Set right hand side Y values
		for (int i = 1; i < n - 1; ++i)
			rhs[i] = 4 * knots[i][1] + 2 * knots[i + 1][1];
		rhs[0] = knots[0][1] + 2 * knots[1][1];
		rhs[n - 1] = (8 * knots[n - 1][1] + knots[n][1]) / 2.0;
		// Get first control points Y-values
		double* y = GetFirstControlPoints(rhs, n);

		// Fill output arrays.
		for (int i = 0; i < n; ++i)
		{
			// First control point
			firstControlPoints[i][0] = x[i];
			firstControlPoints[i][1] = y[i];
			// Second control point
			if (i < n - 1) {
				secondControlPoints[i][0] = 2 * knots[i + 1][0] - x[i + 1];
				secondControlPoints[i][1] = 2 * knots[i + 1][1] - y[i + 1];
			}
			else {
				secondControlPoints[i][0] = (knots[n][0] + x[n - 1]) / 2;
				secondControlPoints[i][1] = (knots[n][1] + y[n - 1]) / 2;
			}
		}
	}

private:
	static double* GetFirstControlPoints(double* rhs, int rhsLength)
	{
		int n = rhsLength;
		double* x = new double[n]; // Solution vector.
		double* tmp = new double[n]; // Temp workspace.

		double b = 2.0;
		x[0] = rhs[0] / b;
		for (int i = 1; i < n; i++) // Decomposition and forward substitution.
		{
			tmp[i] = 1 / b;
			b = (i < n - 1 ? 4.0 : 3.5) - tmp[i];
			x[i] = (rhs[i] - x[i - 1]) / b;
		}
		for (int i = 1; i < n; i++)
			x[n - i - 1] -= tmp[n - i] * x[n - i]; // Backsubstitution.

		return x;
	}
};