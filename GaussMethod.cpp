//
// Created by vadim on 03.03.2021.
//

#include "GaussMethod.h"

#define ZERO_PRECISION 0.0000000000001

static int max_elem_ind(const double * vec, unsigned size) {
    int max_ind = 0;
    for (unsigned i = 0; i < size; ++i) {
        if (vec[i] > vec[max_ind])
            max_ind = static_cast<int>(i);
        if (fabs(0 - vec[i]) < ZERO_PRECISION) {
            max_ind = -1;
            break;
        }
    }
    return max_ind;
}

Matrix GaussMethod::MakeIdentityMatrix(const Matrix &a) {
    Matrix retv(a);

    int max_ind = 0;
    auto tmp_vec = new double[a.sizerow()];

    for (int i = 0; i < a.sizerow(); ++i) {
        retv.GetCol(i, tmp_vec);
        max_ind = max_elem_ind(tmp_vec + i, a.sizerow() - i) + i;
        if (max_ind == -1) {
            break;
        }
        if (max_ind != i)
            retv.rowSwap(max_ind, i);
        retv.multiplyRowByConst(i, 1.0 / retv[i][i]);
        for (unsigned j = 0; j < retv.sizerow(); ++j) {
            if (j == i) continue;
            retv.sumRowByConstAndRow(j, i, -retv[j][i]);
        }
    }
    delete [] tmp_vec;
    return retv;
}