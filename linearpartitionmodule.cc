/*
 * Python interface for LinearPartition
 *
 * Copyright 2024 Hyeshik Chang
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * “Software”), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
 * NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Python.h"
#include "numpy/arrayobject.h"
#include <fstream>
#include <algorithm>
int
trap_fprintf(FILE *fp, const char *fmt, ...)
{
    /* Block outputs to stderr */
    return 0;
}

#define private protected
#undef fprintf
#define fprintf trap_fprintf

// Monkey patch LinearPartition symbols to allow double-linking of E and V

/* EternaFold model */
#define main _linearpartition_e_main
#define BeamCKYParser LPE_BeamCKYParser
#define multi_base LPE_multi_base
#define multi_unpaired LPE_multi_unpaired
#define multi_paired LPE_multi_paired
#define external_unpaired LPE_external_unpaired
#define external_paired LPE_external_paired
#define base_pair LPE_base_pair
#define internal_1x1_nucleotides LPE_internal_1x1_nucleotides
#define helix_stacking LPE_helix_stacking
#define terminal_mismatch LPE_terminal_mismatch
#define bulge_0x1_nucleotides LPE_bulge_0x1_nucleotides
#define helix_closing LPE_helix_closing
#define dangle_left LPE_dangle_left
#define dangle_right LPE_dangle_right
#define internal_explicit LPE_internal_explicit
#define hairpin_length LPE_hairpin_length
#define bulge_length LPE_bulge_length
#define internal_length LPE_internal_length
#define internal_symmetric_length LPE_internal_symmetric_length
#define internal_asymmetry LPE_internal_asymmetry
#define hairpin_length_at_least LPE_hairpin_length_at_least
#define bulge_length_at_least LPE_bulge_length_at_least
#define internal_length_at_least LPE_internal_length_at_least
#define internal_symmetric_length_at_least LPE_internal_symmetric_length_at_least
#define internal_asymmetry_at_least LPE_internal_asymmetry_at_least
#include "contrib/feature_weight_e.h"
#include "LinearPartition.cpp"
#undef multi_base
#undef multi_unpaired
#undef multi_paired
#undef external_unpaired
#undef external_paired
#undef base_pair
#undef internal_1x1_nucleotides
#undef helix_stacking
#undef terminal_mismatch
#undef bulge_0x1_nucleotides
#undef helix_closing
#undef dangle_left
#undef dangle_right
#undef internal_explicit
#undef hairpin_length
#undef bulge_length
#undef internal_length
#undef internal_symmetric_length
#undef internal_asymmetry
#undef hairpin_length_at_least
#undef bulge_length_at_least
#undef internal_length_at_least
#undef internal_symmetric_length_at_least
#undef internal_asymmetry_at_least

#define lpv /* Vienna model */
#undef FASTCKY_BEAMCKYPAR_H
#undef FASTCKY_W
#undef main
#undef BeamCKYParser
#define main _linearpartition_v_main
#define BeamCKYParser LPV_BeamCKYParser
#define hash_pair LPV_hash_pair
#define comp LPV_comp
#define State LPV_State
#define Fast_Exp LPV_Fast_Exp
#define Fast_LogExpPlusOne LPV_Fast_LogExpPlusOne
#define Fast_LogPlusEquals LPV_Fast_LogPlusEquals
#define quickselect LPV_quickselect
#define quickselect_partition LPV_quickselect_partition
#define rtrim LPV_rtrim
#define pf_type LPV_pf_type
#define value_type LPV_value_type
#include "contrib/pseudouridine.h"
#include "contrib/intl11.h"
#include "contrib/intl21.h"
#include "contrib/intl22.h"
#include "LinearPartition.cpp"
#undef State
#undef main
#undef private
#undef fprintf
#define PAIR_TO_NUM(z) \
    ((z) == 5 ? std::make_pair(0, 3) : \
    ((z) == 1 ? std::make_pair(1, 2) : \
    ((z) == 2 ? std::make_pair(2, 1) : \
    ((z) == 3 ? std::make_pair(2, 3) : \
    ((z) == 4 ? std::make_pair(3, 2) : \
    ((z) == 6 ? std::make_pair(3, 0) : \
    std::make_pair(-1, -1)))))))

struct basepair_prob {
    int32_t i;
    int32_t j;
    double prob;
};

static PyArray_Descr *partition_return_descr;

#define GET_BASEPAIR_PROB \
    PyObject * \
    get_basepair_prob(void) \
    { \
        PyArrayObject *res; \
        npy_intp dim; \
        struct basepair_prob *bpp; \
    \
        dim = Pij.size(); \
    \
        res = (PyArrayObject *)PyArray_SimpleNewFromDescr(1, &dim, partition_return_descr); \
        if (res == NULL) \
            return NULL; \
        Py_INCREF(partition_return_descr); \
    \
        assert(partition_return_descr->elsize == sizeof(struct basepair_prob)); \
        bpp = (struct basepair_prob *)PyArray_DATA(res); \
    \
        for (auto it = Pij.begin(); it != Pij.end(); ++it) { \
            bpp->i = it->first.first - 1; \
            bpp->j = it->first.second - 1; \
            bpp->prob = it->second; \
            bpp++; \
        } \
    \
        return (PyObject *)res; \
    }

class EternaBeamCKYParser : public LPE_BeamCKYParser {
public:
    using LPE_BeamCKYParser::LPE_BeamCKYParser;

    GET_BASEPAIR_PROB

    double
    get_free_energy(void)
    {
        State& viterbi=bestC[seq_length - 1];
        return viterbi.alpha;
    }
};

class ViennaBeamCKYParser : public LPV_BeamCKYParser {
public:
    using LPV_BeamCKYParser::LPV_BeamCKYParser;

    GET_BASEPAIR_PROB

    double
    get_free_energy(void)
    {
        LPV_State& viterbi=bestC[seq_length - 1];
        return -kT * viterbi.alpha / 100.0;
    }
};

PyDoc_STRVAR(linearpartition_partition_doc,
"partition(seq)\n\
\n\
Return the base-pairing probability matrix and ensemble free energy \
predicted by LinearPartition.");

// Using modified base pair parameters for pseudouridine
// As we are hot swapping the parameter files we need to keep track of the original values
// By contrast, ViennaRNA uses a soft constraint method to handle modified bases
// TODO: Implement a soft constraint method
bool used_psi = false;
int OriTerminalAU37 = TerminalAU37;

static PyObject *
linearpartition_partition(PyObject *self, PyObject *args, PyObject *kwds)
{
    const char *seq, *engine="vienna", *mod = "none";
    int beamsize=100, dangles=2;
    static const char *kwlist[] = {"seq", "engine", "mod",
                                   "beamsize", "dangles", NULL};
    enum { ETERNA, VIENNA } engine_enum;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|ssii:partition",
                                     (char**)kwlist, &seq, &engine, &mod,
                                     &beamsize, &dangles))
        return NULL;

    if (strcmp(engine, "eterna") == 0)
        engine_enum = ETERNA;
    else if (strcmp(engine, "vienna") == 0)
        engine_enum = VIENNA;
    else {
        PyErr_SetString(PyExc_ValueError,
                        "engine must be either 'eterna' or 'vienna'");
        return NULL;
    }


    if (strcmp(mod, "none") == 0) {
        if (used_psi){
            TerminalAU37 = OriTerminalAU37;
            std::memcpy(stack37, ori_stack37, sizeof(stack37));
            std::memcpy(int11_37, ori_int11_37, sizeof(int11_37));
            std::memcpy(int21_37, ori_int21_37, sizeof(int21_37));
            std::memcpy(int22_37, ori_int22_37, sizeof(int22_37));
            used_psi = false;
        }
    } else if (strcmp(mod, "psi") == 0) {
        if (engine_enum == ETERNA) {
            PyErr_SetString(PyExc_ValueError,
                    "Currently EternaFold does not support modified bases.");
        } else if (engine_enum == VIENNA) {
            TerminalAU37 = ModTerminalAP37;
            if (!used_psi){
                for (int i = 0; i < NBPAIRS + 1; i++) {
                    for (int j = 0; j < NBPAIRS + 1; j++) {
                        stack37[i][j] += diff_psi[i][j];
                        for (int k = 0; k < NOTON; k++){
                            for (int l = 0; l < NOTON; l++){
                                int11_37[i][j][k][l] += diff_psi[i][j];
                                for (int m = 0; m < NOTON; m++){
                                    int21_37[i][j][k][l][m] += diff_psi[i][j];
                                    for (int n = 0; n < NOTON; n++){
                                        int22_37[i][j][k][l][m][n] += diff_psi[i][j];
                                    }
                                }
                            }
                        }
                    }
                }
            used_psi = true;
            }
    } 
    } else {
        PyErr_SetString(PyExc_ValueError,
                    "mod must be 'modified bases'.\n"
                    "Currently only no modifications and 'psi' (pseudouridine) is supported.");
        return NULL;
    }

    string rna_seq(seq);
    PyObject *probmtx;
    double free_energy;

    /* Call LinearPartition */
    switch (engine_enum) {
    case ETERNA: {
        EternaBeamCKYParser parser(beamsize, true, false, "", "", false, 0.0,
            "", false, 3.0, "", false, false, 0.3, "", "", false, dangles);
        Py_BEGIN_ALLOW_THREADS
        parser.parse(rna_seq);
        Py_END_ALLOW_THREADS

        probmtx = parser.get_basepair_prob();
        if (probmtx == NULL)
            return NULL;
        free_energy = parser.get_free_energy();
        break;
    }

    case VIENNA: {
        ViennaBeamCKYParser parser(beamsize, true, false, "", "", false, 0.0,
            "", false, 3.0, "", false, false, 0.3, "", "", false, dangles);
        Py_BEGIN_ALLOW_THREADS
        parser.parse(rna_seq);
        Py_END_ALLOW_THREADS

        probmtx = parser.get_basepair_prob();
        if (probmtx == NULL)
            return NULL;
        free_energy = parser.get_free_energy();
        break;
    }

    default:
        PyErr_SetString(PyExc_RuntimeError, "unknown engine");
        return NULL;
    }

    PyObject *ret=Py_BuildValue("Od", probmtx, free_energy);
    Py_DECREF(probmtx);

    return ret;
}

static PyMethodDef linearpartition_methods[] = {
    {"partition",   (PyCFunction)linearpartition_partition,
     METH_VARARGS | METH_KEYWORDS, linearpartition_partition_doc},
    {NULL,          NULL} /* sentinel */
};

PyDoc_STRVAR(module_doc,
"CPython interface to LinearPartition");

static struct PyModuleDef linearpartitionmodule = {
    PyModuleDef_HEAD_INIT,
    "linearpartition",
    module_doc,
    0,
    linearpartition_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit_linearpartition(void)
{

    import_array();

    PyObject *op=Py_BuildValue("[(s, s), (s, s), (s, s)]",
                               "i", "i4", "j", "i4", "prob", "f8");
    if (op == NULL)
        return NULL;

    if (PyArray_DescrConverter(op, &partition_return_descr) == NPY_FAIL) {
        Py_DECREF(op);
        return NULL;
    }
    Py_DECREF(op);

    return PyModuleDef_Init(&linearpartitionmodule);
}
