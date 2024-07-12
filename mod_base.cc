#include <cstring>
#include <string>
#include <Python.h>
#include "contrib/mod_base.h"
#include "LinearPartition.cpp"
#define PAIR_TO_NUM(z) \
    ((z) == 5 ? std::make_pair(0, 3) : \
    ((z) == 1 ? std::make_pair(1, 2) : \
    ((z) == 2 ? std::make_pair(2, 1) : \
    ((z) == 3 ? std::make_pair(2, 3) : \
    ((z) == 4 ? std::make_pair(3, 2) : \
    ((z) == 6 ? std::make_pair(3, 0) : \
    std::make_pair(-1, -1)))))))
    
using namespace std;

// Placeholder for the original values
int OriTerminalAU37;
int ori_stack37[NBPAIRS + 1][NBPAIRS + 1];
int ori_int11_37[NBPAIRS + 1][NBPAIRS + 1][NOTON][NOTON];
int ori_int21_37[NBPAIRS + 1][NBPAIRS + 1][NOTON][NOTON][NOTON];
int ori_int22_37[NBPAIRS + 1][NBPAIRS + 1][NOTON][NOTON][NOTON][NOTON];
int ori_mismatchI37[NBPAIRS+1][NOTON][NOTON];
int ori_mismatchH37[NBPAIRS+1][NOTON][NOTON];
int ori_mismatchM37[NBPAIRS+1][NOTON][NOTON];
int ori_mismatch23I37[NBPAIRS+1][NOTON][NOTON];
int ori_mismatch1nI37[NBPAIRS+1][NOTON][NOTON];
int ori_mismatchExt37[NBPAIRS+1][NOTON][NOTON];

void update_bases(string mod){

    // Differences with the original values
    int diff[NBPAIRS + 1][NBPAIRS + 1];
    memset(diff, 0, sizeof(diff));
    int ModTerminal = 0;
    
    // Store the original values
    memcpy(ori_stack37, stack37, sizeof(stack37));
    memcpy(ori_int11_37, int11_37, sizeof(int11_37));
    memcpy(ori_int21_37, int21_37, sizeof(int21_37));
    memcpy(ori_int22_37, int22_37, sizeof(int22_37));
    memcpy(ori_mismatchI37, mismatchI37, sizeof(mismatchI37));
    memcpy(ori_mismatchH37, mismatchH37, sizeof(mismatchH37));
    memcpy(ori_mismatchM37, mismatchM37, sizeof(mismatchM37));
    memcpy(ori_mismatch23I37, mismatch23I37, sizeof(mismatch23I37));
    memcpy(ori_mismatch1nI37, mismatch1nI37, sizeof(mismatch1nI37));
    memcpy(ori_mismatchExt37, mismatchExt37, sizeof(mismatchExt37));
    
    if (mod == "psi"){
        ModTerminal = ModTerminalAP37;
        memcpy(diff, diff_psi, sizeof(diff_psi));
    } else if (mod == "none"){
        // No modifications
    } else {
        PyErr_SetString(PyExc_ValueError,
            "mod must be 'modified bases'.\n"
            "Currently only no modifications and 'psi' (pseudouridine) is supported.");
        return;
    }


    if (mod != "none"){
        TerminalAU37 = ModTerminal;
        for (int i = 0; i < NBPAIRS + 1; i++) {
            for (int j = 0; j < NBPAIRS + 1; j++) {
                int n1, n2, m1, m2;
                std::tie(n1, n2) = PAIR_TO_NUM(i);
                std::tie(m1, m2) = PAIR_TO_NUM(j);

                stack37[i][j] += diff[i][j];
                mismatchI37[i][m1][m2] += diff[i][j];
                mismatchI37[j][n1][n2] += diff[i][j];
                mismatchH37[i][m1][m2] += diff[i][j];
                mismatchH37[j][n1][n2] += diff[i][j];
                mismatchM37[i][m1][m2] += diff[i][j];
                mismatchM37[j][n1][n2] += diff[i][j];
                mismatch23I37[i][m1][m2] += diff[i][j];
                mismatch23I37[j][n1][n2] += diff[i][j];
                mismatch1nI37[i][m1][m2] += diff[i][j];
                mismatch1nI37[j][n1][n2] += diff[i][j];
                mismatchExt37[i][m1][m2] += diff[i][j];
                mismatchExt37[j][n1][n2] += diff[i][j];
                
                for (int k = 0; k < NOTON; k++){
                    for (int l = 0; l < NOTON; l++){
                        int11_37[i][j][k][l] += diff[i][j];
                        for (int m = 0; m < NOTON; m++){
                            int21_37[i][j][k][l][m] += diff[i][j];
                            for (int n = 0; n < NOTON; n++){
                                int22_37[i][j][k][l][m][n] += diff[i][j];
                            }
                        }
                    }
                }
            }
        }
    } else if (mod == "none"){
            TerminalAU37 = OriTerminalAU37;
            memcpy(stack37, ori_stack37, sizeof(stack37));
            memcpy(int11_37, ori_int11_37, sizeof(int11_37));
            memcpy(int21_37, ori_int21_37, sizeof(int21_37));
            memcpy(int22_37, ori_int22_37, sizeof(int22_37));
            memcpy(mismatchI37, ori_mismatchI37, sizeof(mismatchI37));
            memcpy(mismatchH37, ori_mismatchH37, sizeof(mismatchH37));
            memcpy(mismatchM37, ori_mismatchM37, sizeof(mismatchM37));
            memcpy(mismatch23I37, ori_mismatch23I37, sizeof(mismatch23I37));
            memcpy(mismatch1nI37, ori_mismatch1nI37, sizeof(mismatch1nI37));
            memcpy(mismatchExt37, ori_mismatchExt37, sizeof(mismatchExt37));
    }
}

