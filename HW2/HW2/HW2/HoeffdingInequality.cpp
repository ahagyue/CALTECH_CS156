#include <iostream>
#include <random>

#include "Matrix.h"

void exp_hoeffding() {

    double nu1 = 0, nu2 = 0, nu3 = 0;

    for (int exp = 0; exp < 100000; exp++) {

        Matrix<int> mat(10, 1000, 0);
        mat.random(0, 1, "uniform");

        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<int> dist(0, 999);

        int c1 = 0;
        int c2 = dist(gen);
        int c3 = 0;

        int min_head = 10;
        for (unsigned i = 0; i < mat.get_col(); i++) {
            int head = 0;
            for (unsigned j = 0; j < mat.get_row(); j++) {
                head += mat(j, i);
            }
            if (head < min_head) {
                c3 = i;
                min_head = head;
            }
        }

        for (unsigned i = 0; i < 10; i++) {
            nu1 += double(mat(i, c1)) / 10;
            nu2 += double(mat(i, c2)) / 10;
            nu3 += double(mat(i, c3)) / 10;
        }

        if (exp % 1000 == 0) {
            cout << nu1 / exp << " " << nu2 / exp << " " << nu3 / exp << endl;
        }
    }

    cout << nu1 / 100000 << " " << nu2 / 100000 << " " << nu3 / 100000 << endl;
}