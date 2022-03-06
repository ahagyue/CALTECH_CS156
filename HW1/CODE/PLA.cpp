// problem 7, 8, 9, 10

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
	int PLA(TargetFunction, Matrix<double>, bool);
};

Hypothesis::Hypothesis() : W(3, 1, 0) {}

Matrix<double> Hypothesis::operator()(const Matrix<double>& p) {
	return ( p * W > 0).typecast<double>() * 2.0 - 1.0;
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
    
    int exp_num = 1000;
    double iterations1 = 0, iterations2 = 0;
    double E1 = 0, E2 = 0;
    bool verbose = false;

    for(int exp = 1; exp <= exp_num; exp++) {
        Matrix<double> x1(10, 2, 0);
        Matrix<double> x2(100, 2, 0);
        Matrix<double> x_out(1000, 2, 0);
        Matrix<double> b1(10, 1, 1);
        Matrix<double> b2(100, 1, 1);
        Matrix<double> b_out(1000, 1, 1);

        x1.random(-1, 1, "uniform");
        x2.random(-1, 1, "uniform");
        x_out.random(-1, 1, "uniform");

        x1.concatenate(b1, 0);
        x2.concatenate(b2, 0);
        x_out.concatenate(b_out, 0);

        TargetFunction target;

        Hypothesis hyp1;
        Hypothesis hyp2;

        iterations1 += (double)(hyp1.PLA(target, x1, false)) / exp_num;
        iterations2 += (double)(hyp2.PLA(target, x2, false)) / exp_num;

        E1 += (1 - accuracy(target, hyp1, x_out)) / exp_num;
        E2 += (1 - accuracy(target, hyp2, x_out)) / exp_num;

        if(exp % 100 == 0 && verbose) {
            cout << "N = 10  : average iteration = " << iterations1 << ", average error rate = " << E1 << endl;
            cout << "N = 100 : average iteration = " << iterations2 << ", average error rate = " << E2 << endl;
        }
    }

    cout << "N = 10  : average iteration = " << iterations1 << ", average error rate = " << E1 << endl;
    cout << "N = 100 : average iteration = " << iterations2 << ", average error rate = " << E2 << endl;
}