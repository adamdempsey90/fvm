#!/usr/bin/env python
from __future__ import print_function
default_pars = {'nx1': int,
        'nx2': int,
        'nx3': int,
        'nscalars': int,
        'gamma' : float,
        'cfl' : float,
        'x1_min': float,
        'x1_max': float,
        'x2_min': float,
        'x2_max': float,
        'x3_min': float,
        'x3_max' : float,
        'tend' : float,
        'Nout' : int,
        'Hout' : int,
        'outputname': str}


def determine_type(val):
    import ast
    if val.lower() == 'no' or val.lower() == 'yes':
        return bool
    try:
        return type(ast.literal_eval(val))
    except ValueError:
        return str

def load_par_file(fname):
    with open(fname,'r') as f:
        lines = [[y.strip() for y in x.split('=')] for x in f.readlines() if x[0]!='#' and len(x.strip())>0 and '=' in x]
    return lines
def read_defs_file(fname):
    with open(fname,'r') as f:
        lines = [x.strip() for x in f.readlines() if '#define' in x and x[0] == '#' and len(x.strip())>0 ]
    return lines


def create_int_block(key,first=False):
    if first:
        out = r"""    if"""
    else:
        out = r"""    else if"""
    out += """ (strcmp(name,"{}")==0) """.format(key)
    out += """ { params->""" + key + """= int_val; PRINT_INT(name,int_val); }
    """

    return out

def create_float_block(key,first=False):
    if first:
        out = r"""    if"""
    else:
        out = r"""    else if"""
    out += """ (strcmp(name,"{}")==0) """.format(key)
    out += """ { params->""" + key + """= double_val; PRINT_DOUBLE(name,double_val); }
    """
    return out

def create_bool_block(key,first=False):
    if first:
        out = r"""    if"""
    else:
        out = r"""    else if"""
    out += """ (strcmp(name,"{}")==0) """.format(key)
    out += """ { params->""" + key + """= bool_val; PRINT_STR(name,str_val); }
    """

    return out
def create_string_block(key,first=False):
    if first:
        out = r"""    if"""
    else:
        out = r"""    else if"""
    out += """ (strcmp(name,"{}")==0) """.format(key)
    out += """ { sprintf(params->""" + key + ""","%s",str_val); PRINT_STR(name,str_val); }
    """
    return out
def create_par_file(lines,fname):

    output = r"""#include "defs.h"
#include <ctype.h>
#define PRINT_DOUBLE(NAME,VAL) printf("\t%s = %lg\n",NAME,VAL)
#define PRINT_INT(NAME,VAL) printf("\t%s = %d\n",NAME,VAL)
#define PRINT_STR(NAME,VAL) printf("\t%s = %s\n",NAME,VAL)
#define FPRINT_DOUBLE(F,NAME,VAL) fprintf(f,"%s = %lg\n",NAME,VAL)
#define FPRINT_INT(F,NAME,VAL) fprintf(f,"%s = %d\n",NAME,VAL)
#define FPRINT_STR(F,NAME,VAL) fprintf(f,"%s = %s\n",NAME,VAL)
void set_var(char *name,int int_val, double double_val, int bool_val, char *str_val, Parameters *params) {
"""

    int_blocks =[]
    bool_blocks=[]
    float_blocks=[]
    str_blocks = []
    for i,line in enumerate(default_pars.items()):
        key,t = line
        key = key.lower()
        kargs = {'first': False}
        if t is int:
            kargs['first'] = len(int_blocks) == 0
            int_blocks.append(create_int_block(key,**kargs))
        elif t is bool:
            bool_blocks.append(create_bool_block(key,**kargs))
        elif t is float:
            float_blocks.append(create_float_block(key,**kargs))
        elif t is str:
            str_blocks.append(create_string_block(key,**kargs))
        else:
            print('{} has no type! Given type '.format(key),t)
    for i,line in enumerate(lines):
        key,val = line
        if key not in default_pars:
            key = key.lower()
            t = determine_type(val)
            kargs = {'first': False}
            if t is int:
                kargs['first'] = len(int_blocks) == 0
                int_blocks.append(create_int_block(key,**kargs))
            elif t is bool:
                bool_blocks.append(create_bool_block(key,**kargs))
            elif t is float:
                float_blocks.append(create_float_block(key,**kargs))
            elif t is str:
                str_blocks.append(create_string_block(key,**kargs))
            else:
                print('{} has no type! Given type '.format(key),t)
    output += '\n'.join(int_blocks)
    output += '\n'.join(float_blocks)
    output += '\n'.join(bool_blocks)
    output += '\n'.join(str_blocks)
    output += '\nreturn;\n}\n'

    output += r"""
void parse_argument(int argc, char *argv[], Parameters *params) {
    int j;
    unsigned int i;
    char name[100],strval[100];
    double dval;
    int ival;
    int bool_val;
    char testbool;


    for(j=0;j<argc;j++) {
        sscanf(argv[j],"%32[^=]=%s",name,strval);
        dval = atof(strval);
        ival = atoi(strval);
        testbool = toupper(strval[0]);
        if (testbool == 'Y') bool_val = TRUE;
        else bool_val = FALSE;
        for (i = 0; i<strlen(name); i++) name[i] = (char)tolower(name[i]);
        set_var(name,ival,dval,bool_val,strval,params);
    }



    return;
}

void read_param_file(char *fname, int argc, char *argv[], Parameters *params) {
    FILE *f;

    char tok[20] = "\t :=>";

    char line[100],name[100],strval[100];
    char *data;
    double temp;
    int status;
    int int_val;
    int bool_val;
    char testbool;
    unsigned int i;

    f= fopen(fname,"r");

    while (fgets(line,100,f)) {
        status = sscanf(line,"%s",name);
        if (name[0] != '#' && status == 1) {

             data = line + (int)strlen(name);
             sscanf(data + strspn(data,tok),"%lf",&temp);
             sscanf(data + strspn(data,tok),"%s",strval);
            int_val = (int)temp;
            testbool = toupper(strval[0]);
            if (testbool == 'Y') bool_val = TRUE;
            else bool_val = FALSE;

            for (i = 0; i<strlen(name); i++) name[i] = (char)tolower(name[i]);

            set_var(name,int_val,temp,bool_val,strval,params);

        }
    }

    if (argc > 0) {
        printf("Redefined on the command line:\n");
        parse_argument(argc,argv,params);
    }

    return;
}


    """

    with open(fname,'w') as f:
        f.write(output)

def create_struct(lines):
    out_lines = ['#ifdef ISFLOAT\n#define real float\n#else\n#define real double\n#endif\n',
            'typedef struct Parameters {']
    for key,t in default_pars.items():
        key = key.lower()
        if t is int or t is bool:
            out_lines.append('\tint {};'.format(key))
        elif t is float:
            out_lines.append('\treal {};'.format(key))
        elif t is str:
            out_lines.append('\tchar {}[512];'.format(key))
        else:
            print('{} has no type! Given type '.format(key),t)

    for line in lines:
        key,val = line
        if key not in default_pars:
            key = key.lower()
            t = determine_type(val)
            if t is int or t is bool:
                out_lines.append('\tint {};'.format(key))
            elif t is float:
                out_lines.append('\treal {};'.format(key))
            elif t is str:
                out_lines.append('\tchar {}[512];'.format(key))
            else:
                print('{} has no type! Given type '.format(key),t)

    out_lines.append('} Parameters;\n');


    out_lines.append('void read_param_file(char *fname, int argc, char *argv[], Parameters *params);\n')

    return out_lines




if __name__ == "__main__":
    """
        Pass the directory containing the initialization file,
        parameter file, and problem defs.

        Example for implosion test located in

          src/tests/2D/imp/
          src/tests/2D/imp/imp.cu  --> init file
          src/tests/2D/imp/imp.par --> parameters
          src/tests/2D/imp/imp.h   --> problem defs

        Configure and compile this test using,

        ./configure -prob src/tests/2D/imp
        make
    """
    import argparse
    import shutil
    parser = argparse.ArgumentParser()
    parser.add_argument('-prob',type=str,default='src/tests/2D/imp',help='Problem directory')
    parser.add_argument('--prof',action='store_true',help='Enabling profiling. No outputs will be written.')
    parser.add_argument('--silent',action='store_true',help='Silences all output to stdout.')


    args = vars(parser.parse_args())

    directory = args['prob']
    prof = args['prof']
    silent = args['silent']

    if directory[-1] != '/':
        directory += '/'

    problem_name = directory.split('/')[-2]

    parfile = directory + problem_name + '.par'
    defsfile = directory + problem_name + '.h'
    initfile = directory + problem_name + '.cu'



    shutil.copy(defsfile,'src/prob.h')
    shutil.copy(initfile,'src/prob.cu')





    lines = load_par_file(parfile)
    create_par_file(lines,'src/read_pars.c')


    struct_lines = create_struct(lines)
    defs_lines = read_defs_file(defsfile)
    if prof:
        defs_lines.append("#define PROF")
    if silent:
        defs_lines.append("#define SILENT")

    with open('src/prob.h','w') as f:

        f.write('\n'.join(defs_lines) + '\n\n')
        f.write('\n'.join(struct_lines) + '\n')
