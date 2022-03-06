//problem  8, 9, 10

#include "Matrix.h"

// transformation function

Matrix<double> secondary_nonlinearity(Matrix<double> points) {
    Matrix<double> transformed(points.get_row(), points.get_col() * (points.get_col() + 1) / 2, 0);
    
    for(unsigned i = 0; i < points.get_row(); i++) {
        int l = 0;
        for(unsigned j = 0; j < points.get_col(); j++) {
            for(unsigned k = j; k < points.get_col(); k++) {
                transformed(i, l) = points(i, j) * points(i, k);
                l++;
            }
        }
    }
    return transformed;
}

//target function
Matrix<double> TargetFunction(Matrix<double> points) {
    points = secondary_nonlinearity(points);
    Matrix<double> W(6, 1, 0);
    W(0, 0) = 1;
    W(3, 0) = 1;
    W(5, 0) = -0.6;

    return (points * W > 0.0) * 2.0 - 1.0;
}

//-------------hypothesis------------
class Hypothesis {
private:
	Matrix<double> W;
    bool nonlinear;
    Matrix<double> (*nonlinear_transformation)(Matrix<double>);
public:
	Hypothesis();
    Hypothesis(int);
    Hypothesis(int, Matrix<double> (*nonlinear_transformation)(Matrix<double>));

	Matrix<double> operator()(const Matrix<double>&);
    Matrix<double> get_W();

	void linear_regression(const Matrix<double>&, const Matrix<double>&);
};

Hypothesis::Hypothesis() : W(3, 1, 0) {
    nonlinear = false;
	W.random(-1, 1, "uniform");
}

Hypothesis::Hypothesis(int w_size) : W(w_size, 1, 0) {
    nonlinear = false;
    W.random(-1, 1, "uniform");
}

Hypothesis::Hypothesis(int w_size, Matrix<double> (*nonlinear_transformation)(Matrix<double>)) : W(w_size, 1, 0) {
    nonlinear = true;
	W.random(-1, 1, "uniform");
    this -> nonlinear_transformation = nonlinear_transformation;
}

Matrix<double> Hypothesis::operator()(const Matrix<double>& p) {
    if(nonlinear) {
	    return ( nonlinear_transformation(p) * W > 0).typecast<double>() * 2.0 - 1.0;
    }
	return ( p * W > 0).typecast<double>() * 2.0 - 1.0;
}

Matrix<double> Hypothesis::get_W() {
    return W;
}

void Hypothesis::linear_regression(const Matrix<double>& x, const Matrix<double>& y) {
    if(nonlinear) {
        Matrix<double> nonlinear_x = nonlinear_transformation(x);
        W = (nonlinear_x.inverse() * y);
    } else {
        W = (x.inverse() * y);
    }
}
// --------- hypothesis----------


// accuracy function
double accuracy(Hypothesis hyp, Matrix<double> x, Matrix<double> target_y) {
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

int main() {
    
    int exp_num = 1000;
    double linear_E_in = 0;
    double nonlinear_E_in = 0;
    double nonlinear_E_out = 0;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dist(0, 999);

    for(int exp = 1; exp <= exp_num; exp++) {
        // training set generation
        Matrix<double> x(1000, 2, 0); // (x1, x2)
        Matrix<double> b(1000, 1, 1); // bias
        x.random(-1, 1, "uniform");
        x.concatenate(b, 0);
        
        // generation of noisy 'y'
        Matrix<double> y = TargetFunction(x);

        // if you want probabilistic correction due to ignoring overlapped selection, "y.get_row() / 10" -->  "y.get_row() / 9"
        for(int i = 0; i < y.get_row() / 10; i++) {
            int selected = dist(gen);
            y(selected, 0) = 1 - y(selected, 0);    // ignore the case of overlapped selection
        }

        // training set : x, y  | target function : TargetFunction() | hypothesis : Hypothesis class

        // linear hypothesis : problem 8
        Hypothesis linear_hyp;
        linear_hyp.linear_regression(x, y);
        linear_E_in += (1 - accuracy(linear_hyp, x, y)) / exp_num;

        //non-linear hypothesis : problem 9
        Hypothesis nonlinear_hyp(6, &secondary_nonlinearity);
        nonlinear_hyp.linear_regression(x, y);
        nonlinear_E_in += (1 - accuracy(nonlinear_hyp, x, y)) / exp_num;

        if(exp == 1)
            cout << "non-linear regression result (x1^2, x1x2, x1, x2^2, x2, 1) : " << nonlinear_hyp.get_W().T() << endl;

        // problem 10
        // testing set generation
        Matrix<double> x_out(1000, 2, 0); // (x1, x2)
        Matrix<double> b_out(1000, 1, 1); // bias
        x_out.random(-1, 1, "uniform");
        x_out.concatenate(b_out, 0);
        
        // generation of noisy 'y_out'
        Matrix<double> y_out = TargetFunction(x_out);
        for(int i = 0; i < y_out.get_row() / 10; i++) {
            int selected = dist(gen);
            y_out(selected, 0) = 1 - y_out(selected, 0);    // ignore the case of overlapped selection
        }

        nonlinear_E_out += (1 - accuracy(nonlinear_hyp, x_out, y_out)) / exp_num;
    }
    cout << "linear hypothesis E_in      : " << linear_E_in << endl;
    cout << "non-linear hypothesis E_in  : " << nonlinear_E_in << endl;
    cout << "non-linear hypothesis E_out : " << nonlinear_E_out << endl;

}