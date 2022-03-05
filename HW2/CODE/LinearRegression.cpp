//problem 5, 6, 7

#include "Matrix.h"

// --------- target function (linear) ----------
class TargetFunction {
private:
	Matrix<double> W;
public:
	TargetFunction(); // initialize line with random two points in [-1, 1] X [-1, 1] space
	Matrix<double> operator()(const Matrix<double>&);
};

TargetFunction::TargetFunction() : W(3, 1, 0){
	Matrix<double> p1(1, 2, 0);
	Matrix<double> p2(1, 2, 0);

	p1.random(-1., 1., "uniform");
	p2.random(-1., 1., "uniform");

	W(0, 0) = p2(0, 1) - p1(0, 1);
	W(1, 0) = p2(0, 0) - p1(0, 0);
	W(2, 0) = p2(0, 0) * p1(0, 1) - p1(0, 0) * p2(0, 1);
}

Matrix<double> TargetFunction::operator()(const Matrix<double>& p) {
	return ( p * W > 0).typecast<double>() * 2.0 - 1.0;
}

// --------- target function ----------


// --------- hypothesis----------
class Hypothesis {
private:
	Matrix<double> W;
public:
	Hypothesis();
	Matrix<double> operator()(const Matrix<double>&);
	void linear_regression(const Matrix<double>&, const Matrix<double>&);
	int PLA(TargetFunction, Matrix<double>, bool);
};

Hypothesis::Hypothesis() : W(3, 1, 0) {
	W.random(-1, 1, "uniform");
}

Matrix<double> Hypothesis::operator()(const Matrix<double>& p) {
	return ( p * W > 0).typecast<double>() * 2.0 - 1.0;
}

void Hypothesis::linear_regression(const Matrix<double>& x, const Matrix<double>& y) {
	W = (x.inverse() * y);
}

// accuracy function
double accuracy(TargetFunction target, Hypothesis hyp, Matrix<double> x) {
	Matrix<double> target_y = target(x);
	Matrix<double> hyp_y = hyp(x);

	double accuracy = 0;

	for(unsigned i = 0; i < target_y.get_row(); i++) {
		if(target_y(i, 0) == hyp_y(i, 0))
			accuracy += 1;
	}
	accuracy /= target_y.get_row();

	return accuracy;
}

// --------------------

int Hypothesis::PLA(TargetFunction target, Matrix<double> x, bool verbose) {
	int iterations = 0;

	while(accuracy(target, (*this), x) != 1 ) {
		iterations++;
		for(unsigned i = 0; i < x.get_row(); i++) {
			Matrix<double> point = x.get_line(i, 0);
			if(target(point)(0, 0)  != (*this)(point)(0, 0))
				W += target(point)(0, 0) * point.T();
		}

		if(verbose) {
			cout << "iterations : " << iterations << " | accuracy : " << accuracy(target, (*this), x) << " | W : " << W.T() << endl;
		}
	} 

	return iterations;
}
// --------- hypothesis----------

int main() {

	//problem 5, 6

	double E_in = 0, E_out = 0;
	int exp_num = 1000;

	for(int exp = 1; exp <= exp_num; exp++ ) {
		Matrix<double> x_in(100, 3, 0);
		Matrix<double> x_out(1000, 3, 0);
		x_in.random(-1, 1, "uniform");
		x_out.random(-1, 1, "uniform");

		for(unsigned i = 0; i < x_in.get_row(); i++)
			x_in(i, 2) = 1;
		for(unsigned i = 0; i < x_out.get_row(); i++)
			x_out(i, 2) = 1;
			
		TargetFunction target;
		Matrix<double> y = target(x_in);

		Hypothesis hyp;
		hyp.linear_regression(x_in, y);

		E_in += (1 - accuracy(target, hyp, x_in)) / exp_num;
		E_out += (1 - accuracy(target, hyp, x_out)) / exp_num;
	}

	cout << "E_in :  " << E_in << endl;
	cout << "E_out : " << E_out << endl;


	// problem 7
	exp_num = 1000;
	double iterations = 0;
	for(int exp = 1; exp <= exp_num; exp++) {
		Matrix<double> x(10, 3, 0);
		x.random(-1, 1, "uniform");
		for(unsigned i = 0; i < x.get_row(); i++)
			x(i, 2) = 1;
		
		TargetFunction target;
		Matrix<double> y = target(x);

		Hypothesis hyp;
		hyp.linear_regression(x, y);

		iterations += (double)(hyp.PLA(target, x, false)) / exp_num;
	}
	cout << "PLA iterations : " << iterations <<endl;
}