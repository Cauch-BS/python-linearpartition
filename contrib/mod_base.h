/* Original json file for modified base parameters from ViennaRNA
Authors : Graham A. Hudson, Richard J. Bloomingdale, and Brent M. Znosko
Title : Thermodynamic contribution and nearest-neighbor parameters of pseudouridine-adenosine base pairs in oligoribonucleotides
Journal : RNA 19:1474-1482
Year : 2013
DOI : 10.1261/rna.039610.113
{
    "name" : "Pseudouridine",
    "unmodified" : "U",
    "pairing_partners" : [
        "A"
    ],
    "one_letter_code" : "P",
    "fallback" : "U",
    "stacking_energies" : {   --> cbsheen: modifed values on lines 109-120
        "APUA" :  -2.8,       --> cbsheen: "AP" & "UA"
        "CPGA" : -2.77,       --> cbsheen: "AP" & "GC"
        "GPCA" : -3.29,       --> cbsheen: "AP" & "CG"
        "UPAA" : -1.62,       --> cbsheen: "AP" & "AU"
        "PAAU" : -2.10,       --> cbsheen: "PA" & "AU"
        "PCAG" : -2.49,       --> cbsheen: "PA" & "CG"
        "PGAC" : -2.2,        --> cbsheen: "PA" & "GC"
        "PUAA" : -2.74        --> cbsheen: "PA" & "UA"
    },
    "stacking_enthalpies" : {
        "APUA" : -22.08,
        "CPGA" : -16.23,
        "GPCA" : -24.07,
        "UPAA" : -20.81,
        "PAAU" : -12.47,
        "PCAG" : -17.29,
        "PGAC" : -11.19,
        "PUAA" : -26.94
    },
    "terminal_energies" : {   --> cbsheen: modified value on line 70
        "PA" : 0.31,
        "AP" : 0.31
    },
    "terminal_enthalpies" : {
        "PA" : -2.04,
        "AP" : -2.04
    },
}
*/

int ModTerminalAP37 = 31;

int diff_psi[NBPAIRS+1][NBPAIRS+1] = 
    {
    //                             CG          GC          GU          UG          AP          PA          NN   
            {           0,          0,          0,          0,          0,          0,          0,          0}
     /*CG*/,{           0,          0,          0,          0,          0,        -67,        -10,          0}
     /*GC*/,{           0,          0,          0,          0,          0,       -109,        -67,          0}
     /*GU*/,{           0,          0,          0,          0,          0,          0,          0,          0}
     /*UG*/,{           0,          0,          0,          0,          0,          0,          0,          0}
     /*AP*/,{           0,        -67,       -109,          0,          0,       -340,       -256,          0}
     /*PA*/,{           0,        -10,         -9,          0,          0,       -256,       -160,          0}
     /*NN*/,{           0,          0,          0,          0,          0,          0,          0,          0}
    } ;

