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
int
trap_fprintf(FILE *fp, const char *fmt, ...)
{
    /* Block outputs to stderr */
    return 0;
}

#define main _linearpartition_main
#define lpv /* Vienna model */
#define private protected
#undef fprintf
#define fprintf trap_fprintf
#include "LinearPartition.cpp"
#undef main
#undef private
#undef fprintf

struct basepair_prob {
    int32_t i;
    int32_t j;
    double prob;
};

static PyArray_Descr *partition_return_descr;

class MyBeamCKYParser : public BeamCKYParser {
public:
    MyBeamCKYParser(int beamsize, bool no_sharpturn, bool verbose,
                    const string &bpp_file, const string &bpp_file_index,
                    bool pf_only, float bpp_cutoff, const string &forest_file,
                    bool mea, float MEA_gamma, const string &MEA_file_index,
                    bool MEA_bpseq, bool ThreshKnot, float ThreshKnot_threshold,
                    const string &ThreshKnot_file_index, const string &shape_file_path,
                    bool fasta, int dangles)
        : BeamCKYParser(beamsize, no_sharpturn, verbose, bpp_file, bpp_file_index,
                        pf_only, bpp_cutoff, forest_file, mea, MEA_gamma, MEA_file_index,
                        MEA_bpseq, ThreshKnot, ThreshKnot_threshold, ThreshKnot_file_index,
                        shape_file_path, fasta, dangles)
    {
    }

    PyObject *
    get_basepair_prob(void)
    {
        PyArrayObject *res;
        npy_intp dim;
        struct basepair_prob *bpp;

        dim = Pij.size();

        res = (PyArrayObject *)PyArray_SimpleNewFromDescr(1, &dim, partition_return_descr);
        if (res == NULL)
            return NULL;
        Py_INCREF(partition_return_descr);

        assert(partition_return_descr->elsize == sizeof(struct basepair_prob));
        bpp = (struct basepair_prob *)PyArray_DATA(res);

        for (auto it = Pij.begin(); it != Pij.end(); ++it) {
            bpp->i = it->first.first - 1;
            bpp->j = it->first.second - 1;
            bpp->prob = it->second;
            bpp++;
        }

        return (PyObject *)res;
    }

    double
    get_free_energy(void)
    {
        State& viterbi=bestC[seq_length - 1];
        return -kT * viterbi.alpha / 100.0;
    }
};

PyDoc_STRVAR(linearpartition_partition_doc,
"partition(seq)\n\
\n\
Return the base-pairing probability matrix and ensembl free energy \
predicted by LinearPartition.");

static PyObject *
linearpartition_partition(PyObject *self, PyObject *args)
{
    const char *seq;
    Py_ssize_t len;

    if (!PyArg_ParseTuple(args, "s#:partition", &seq, &len))
        return NULL;

    /* LinearPartition arguments */
    int beamsize = 100;
    bool sharpturn = false;
    bool pf_only = false;
    float bpp_cutoff = 0.0;

    float MEA_gamma = 3.0;
    float ThreshKnot_threshold = 0.3;
    int dangles = 2;
    // End of LinearPartion parameters

    string rna_seq(seq);

    /* Call LinearPartition */
    MyBeamCKYParser parser(beamsize, !sharpturn, false, "", "",
                           pf_only, bpp_cutoff, "", false, MEA_gamma, "",
                           false, false, ThreshKnot_threshold, "",
                           "", false, dangles);
    Py_BEGIN_ALLOW_THREADS
    parser.parse(rna_seq);
    Py_END_ALLOW_THREADS

    PyObject *ret, *probmtx;

    probmtx = parser.get_basepair_prob();
    if (probmtx == NULL)
        return NULL;

    ret = Py_BuildValue("Od", probmtx, parser.get_free_energy());
    Py_DECREF(probmtx);

    return ret;
}

static PyMethodDef linearpartition_methods[] = {
    {"partition",   linearpartition_partition,  METH_VARARGS,
     linearpartition_partition_doc},
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
