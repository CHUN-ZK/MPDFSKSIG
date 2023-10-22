#include "espresso.h"
#include "strdup.h"
#include <string.h>
#include <stdio.h>

void getPLA(int argc, char** argv, pPLA* PLA)
{
    int input_type = FD_type; //default
    int out_type = F_type;

    FILE* fp;
    int needs_dcset, needs_offset;
    char* fname;
    int planum = 1;

    //if (argc > 2) {
    //    fp = stdin;
    //    fname = "(stdin)";
    //}
    //else {
        fname = argv[planum];
        if ((fp = fopen(argv[planum], "r")) == NULL) {
            fprintf(stderr, "%s: Unable to open %s\n", argv[0], fname);
            exit(1);
        }
    //}

    needs_dcset = TRUE; //default 
    needs_offset = TRUE;

    if (read_pla(fp, needs_dcset, needs_offset, input_type, PLA) == EOF) {
        fprintf(stderr, "%s: Unable to find PLA on file %s\n", argv[0], fname);
        exit(1);
    }
    (*PLA)->filename = _strdup(fname);
    filename = (*PLA)->filename;
}
