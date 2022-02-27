#include "Matrix.h"

class TargetFunction {
private:
	Matrix<double> W;
public:
	TargetFunction();
	Matrix<double> operator()(const Matrix<double>&);
};

TargetFunction::TargetFunction() {
	Matrix<double> p1(1, 2, 0);
	Matrix<double> p2(1, 2, 0);

	p1.random(-1., 1., "uniform");
	p2.random(-1., 1., "uniform");

	Matrix<double> w(1, 3, 0);
	w(0, 0) = p2(0, 1) - p1(0, 1);
	w(0, 1) = p2(0, 0) - p1(0, 0);
	w(0, 2) = p2(0, 0) * p1(0, 1) - p1(0, 0) * p2(0, 1);
	W = w;
}

Matrix<double> TargetFunction::operator()(const Matrix<double>& p) {
	return W * p;
}

class Hypothesis {
private:
	Matrix<double> W;
public:
	Hypothesis();
	Matrix<double> operator()(const Matrix<double>&);
	void linear_regression();
};

Hypothesis::Hypothesis()
{
}

Matrix<double> Hypothesis::operator()(const Matrix<double>&)
{
	return Matrix<double>();
}

void Hypothesis::linear_regression()
{
}
