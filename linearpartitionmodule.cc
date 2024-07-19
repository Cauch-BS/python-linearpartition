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
#include "ATen/ATen.h"
#include "torch/python.h"
#include "torch/extension.h"
#include <fstream>
// FLAG
// for debuging purposes
// #include <sstream>
// #include <string>

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
#include "contrib/mod_base.h"
#include "LinearPartition.cpp"
#undef State
#undef main
#undef private
#undef fprintf

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

// Using modified base pair parameters for pseudouridine
// As we are hot swapping the parameter files we need to keep track of the original values
// By contrast, ViennaRNA uses a soft constraint method to handle modified bases
// TODO: Implement a soft constraint method

bool used_psi = false; /*True if mod is used. If false, update_bases is not run*/
bool used_m1psi = false; /*True if mod is used. If false, update_bases is not run*/
bool updated_stack = false; /*True if stack37 has been updated*/

double OriTerminalAU37 = TerminalAU37;
double ori_stack37[NBPAIRS + 1][NBPAIRS + 1];
size_t size_of_stack = (NBPAIRS + 1) * (NBPAIRS + 1) * sizeof(double);


/* NOTICE: update_bases should only be run once per modification
Else the return to origin will not work properly.
This is managed with the used_psi and used_m1psi booleans*/
void update_bases(string mod){

    // Differences with the original values
    double diff[NBPAIRS + 1][NBPAIRS + 1];
    memset(diff, 0, size_of_stack);
    double ModTerminal = 0;
    
    // Load the values for modified bases
    if (mod == "psi"){
        ModTerminal = ModTerminalAP37;
        memcpy(diff, diff_psi, size_of_stack);
    } else if (mod == "m1psi"){
        ModTerminal = ModTerminalAP37;
        memcpy(diff, diff_m1psi, size_of_stack);
    }  else if (mod == "none"){
        // No modifications
    } else {
        PyErr_SetString(PyExc_ValueError,
            "mod must be 'modified bases'.\n"
            "Currently only no modifications and 'psi' (pseudouridine) and \n"
            "'m1psi' (N1-methylpseudouridine) is supported.");
        return;
    }

    if (mod != "none"){
        // Store the original values
        if (!used_psi && !used_m1psi && !updated_stack){
            memcpy(ori_stack37, stack37, size_of_stack);
        } else {
            memcpy(stack37, ori_stack37, size_of_stack);
        }
        // Update the values
        TerminalAU37 = ModTerminal;
        for (int i = 0; i < NBPAIRS + 1; i++) {
            for (int j = 0; j < NBPAIRS + 1; j++) {
                stack37[i][j] += diff[i][j];
            }
        }
    } else if (mod == "none"){
            TerminalAU37 = OriTerminalAU37;
            // Return to the original values
            memcpy(stack37, ori_stack37, size_of_stack);
    }
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

static PyObject *
linearpartition_partition(PyObject *self, PyObject *args, PyObject *kwds)
{
    const char *seq, *engine="vienna", *mod = "none";
    int beamsize=100, dangles=2;
    double update_terminal = OriTerminalAU37;
    PyObject *update_stack = NULL;
    static const char *kwlist[] = {"seq", "engine", "mod",
                                   "beamsize", "dangles", "update_terminal",
                                   "update_stack", NULL};
    enum { ETERNA, VIENNA } engine_enum;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|ssiidO:partition",
                                     (char**)kwlist, &seq, &engine, &mod,
                                     &beamsize, &dangles,
                                     &update_terminal, &update_stack))
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

    if (update_stack != NULL || update_terminal != OriTerminalAU37){ 
        if (THPVariable_Check(update_stack)) {
            if ((strcmp(mod, "custom") != 0)&&(strcmp(mod, "none") != 0)){
                PyErr_SetString(PyExc_ValueError,
                            "mod must be 'custom' when update is passed.");
                return NULL;
            }
            mod = "custom";
        } else {
            PyErr_SetString(PyExc_ValueError,
                            "update_stack must be a Pytorch Tensor");
            return NULL;
        }
    }


    if (strcmp(mod, "none") == 0) {
        if ((used_psi || used_m1psi || updated_stack)){
            update_bases("none");
            used_psi = false;
            used_m1psi = false;
            updated_stack = false;
        }
    } else if (strcmp(mod, "psi") == 0) {
        if (engine_enum == ETERNA) {
            PyErr_SetString(PyExc_ValueError,
                    "Currently EternaFold does not support modified bases.");
        } else if (engine_enum == VIENNA) {
            if (!used_psi){
                update_bases("psi");
                used_psi = true;
                used_m1psi = false;
                updated_stack = false;
            }
        } 
    } else if (strcmp(mod, "m1psi") == 0) {
        if (engine_enum == ETERNA) {
            PyErr_SetString(PyExc_ValueError,
                    "Currently EternaFold does not support modified bases.");
        } else if (engine_enum == VIENNA) {
            if (!used_m1psi){
                update_bases("m1psi");
                used_m1psi = true;
                used_psi = false;
                updated_stack = false;
            }
        } 

    } else if (strcmp(mod, "custom") == 0) {
        if (engine_enum == ETERNA) {
            PyErr_SetString(PyExc_ValueError,
                    "Currently EternaFold does not support modified bases.");
        } else if (engine_enum == VIENNA) {
            if ((update_stack != NULL) &&  (THPVariable_Check(update_stack))){
                if (!updated_stack && !used_psi && !used_m1psi){
                    memcpy(ori_stack37, stack37, size_of_stack);
                } else {
                    memcpy(stack37, ori_stack37, size_of_stack);
                }

                auto tensor = THPVariable_Unpack(update_stack);
                if (tensor.dim() != 2){
                    PyErr_SetString(PyExc_ValueError,
                                "update_stack must be a 2D tensor.");
                    return NULL;
                }
                if (tensor.size(0) != NBPAIRS + 1 || tensor.size(1) != NBPAIRS + 1){
                    PyErr_SetString(PyExc_ValueError,
                                "update_stack is of incorrect shape. Reshape to 8 x 8 tensor.");
                    return NULL;
                }
                if (tensor.scalar_type() != torch::kDouble){
                    PyErr_SetString(PyExc_ValueError,
                                "update_stack must be of dtype torch.float64 (double).");
                    return NULL;
                }

                tensor = tensor.contiguous();
                auto accessor = tensor.accessor<double, 2>();
                for (int i = 0; i < NBPAIRS + 1; i++) {
                    for (int j = 0; j < NBPAIRS + 1; j++) {
                        stack37[i][j] += accessor[i][j];
                    }
                }
            } else {
                PyErr_SetString(PyExc_ValueError,
                        "when mod is 'custom', update_stack must be passed as a Pytorch Tensor");
                return NULL;
            }

            TerminalAU37 = update_terminal;

            updated_stack = true;
            used_m1psi = false;
            used_psi = false;
            
        }
    } else {
        PyErr_SetString(PyExc_ValueError,
                        "mod must be 'modified bases'.\n"
                        "Currently only no modifications and 'psi' (pseudouridine) and \n"
                        "'m1psi' (N1-methylpseudouridine) is supported.");
        return NULL;
    }

    // if ( true ){
    //     // For debugging
    //     // Print the Updated Stack
    //     std::ostringstream oss_stack, oss_terminal;
    //     oss_stack << "Stack array contents:\n" << std::setprecision(6) << std::fixed;
    //     for (int i = 0; i < NBPAIRS + 1; i++) {
    //         for (int j = 0; j < NBPAIRS + 1; j++) {
    //             oss_stack << stack37[i][j] << " ";
    //         }
    //         oss_stack << "\n";
    //     }
        
    //     oss_terminal << "Terminal contents:\n" << TerminalAU37 << "\n";
        
    //     std::cerr << oss_stack.str() << std::endl;
    //     std::cerr << oss_terminal.str() << std::endl;
    // }

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
