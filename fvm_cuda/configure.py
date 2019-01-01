#!/usr/bin/env python
from __future__ import print_function
default_pars = {'nx1': (int,128),
        'nx2': (int,1),
        'nx3': (int,1),
        'nscalars': (int,0),
        'gamma' : (float,1.4),
        'cfl' : (float,.2),
        'x1_min': (float,0.),
        'x1_max': (float,1.),
        'x2_min': (float,0.),
        'x2_max': (float,1.),
        'x3_min': (float,0.),
        'x3_max' : (float,1.),
        'tend' : (float,1.),
        'nout0d' : (int,-1),
        'nout1d' : (int,-1),
        'nout2d' : (int,-1),
        'nout3d' : (int,-1),
        'maxsteps' : (int,"1e64"),
        'outputdir' : (str,'out/'),
        'outputname': (str,'fld')}


def default_lines():

    lines = "void default_pars(Parameters *params) {\n"
    for key,val in default_pars.items():
        t,default = val
        if key == 'maxsteps':
            lines += "    params->{} = (long long int){};\n".format(key,str(default))
        elif t != str:
            lines += "    params->{} = {};\n".format(key,str(default))
        else:
            lines += '    strcpy(params->{},"{}");\n'.format(key,default)

    lines += '\n}\n'
    return lines

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
    lines = [ [x[0].lower(), x[1]] for x in lines ]
    return lines
def read_defs_file(fname):
    with open(fname,'r') as f:
        lines = [x.strip() for x in f.readlines() if '#define' in x and x[0] == '#' and len(x.strip())>0 ]
    return lines


def create_int_block(key,first=False):
    if key == 'maxsteps':
        if first:
            out = r"""    if"""
        else:
            out = r"""    else if"""
        out += """ (strcmp(name,"{}")==0) """.format(key)
        out += """ { params->""" + key + """= (long long int)double_val; PRINT_DOUBLE(name,double_val); }
        """
        return out

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
        key,val = line
        t,default = val
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
        key = key.lower()
        print(line, key, val)
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

    output += default_lines()

    with open(fname,'w') as f:
        f.write(output)

def create_struct(lines):
    out_lines = ['#ifdef ISFLOAT\n#define real float\n#else\n#define real double\n#endif\n',
            'typedef struct Parameters {']
    for key,val in default_pars.items():
        key = key.lower()
        t,default = val
        if t is int or t is bool:
            if key == 'maxsteps':
                out_lines.append('\tlong long int {};'.format(key))
            else:
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
    out_lines.append('void default_pars(Parameters *params);\n')

    return out_lines



def add_extra_def(name,args,defs_lines,extra_defs,defname=None):
    if args[name.lower()]:
        if not any([name.lower() in x.lower() for x in defs_lines]):
            extra_defs.append('#define {}'.format(name.upper() if defname is None else defname))
    return

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
    from subprocess import call
    parser = argparse.ArgumentParser()
    parser.add_argument('-prob',type=str,default='src/tests/2D/imp',help='Problem directory')
    parser.add_argument('--prof',action='store_true',help='Enabling profiling. No outputs will be written.')
    parser.add_argument('--silent',action='store_true',help='Silences all output to stdout.')
    parser.add_argument('--pcm',action='store_true',help='Piecewise constant reconstruction.')
    parser.add_argument('--plm',action='store_true',help='Piecewise linear reconstruction.')
    parser.add_argument('--ppm',action='store_true',help='Piecewise parabolic reconstruction.')
    parser.add_argument('--ctu',action='store_true',help='Use CTU algorithm.')
    parser.add_argument('--conduction',action='store_true',help='Enable heat conduction.')
    parser.add_argument('--viscosity',action='store_true',help='Enable viscosity.')
    parser.add_argument('--potential',action='store_true',help='Enable static gravitational potential.')
    parser.add_argument('--hll',action='store_true',help='Use HLL Riemann solver.')
    parser.add_argument('--hllc',action='store_true',help='Use HLLC Riemann solver.')
    parser.add_argument('--exact',action='store_true',help='Use exact Riemann solver.')
    parser.add_argument('--float',action='store_true',help='Use floats instead of doubles.')
    parser.add_argument('--dims1',action='store_true',help='1D problem.')
    parser.add_argument('--dims2',action='store_true',help='2D problem.')
    parser.add_argument('--dims3',action='store_true',help='3D problem.')




    args = vars(parser.parse_args())


    directory = args['prob']




    if directory[-1] != '/':
        directory += '/'

    problem_name = directory.split('/')[-2]

    parfile = directory + problem_name + '.par'
    defsfile = directory + problem_name + '.h'
    initfile = directory + problem_name + '.cu'








    defs_lines = read_defs_file(defsfile)
    extra_defs = []

    # Number of dims
    dims1 = int(args['dims1'])
    dims2 = int(args['dims2'])
    dims3 = int(args['dims3'])
    rsum = dims1 + dims2 + dims3
    if rsum > 0:
        if rsum > 1:
            print('Can only have one of dims1, dims2, dims3 defined!')
            exit()
        else:
            # Check for already defined Riemann solver
            defs_lines = list(filter(lambda x: not any([c in x.lower() for c in ['dims1','dims2','dims3']]), defs_lines))
            if dims1:
                extra_defs.append('#define DIMS1')
            if dims2:
                extra_defs.append('#define DIMS2')
            if dims3:
                extra_defs.append('#define DIMS3')

    # Riemann solver
    hll = int(args['hll'])
    hllc = int(args['hllc'])
    exact = int(args['exact'])
    rsum = hll + hllc + exact
    if rsum > 0:
        if rsum > 1:
            print('Can only have one of HLL, HLLC, EXACT defined!')
            exit()
        else:
            # Check for already defined Riemann solver
            defs_lines = list(filter(lambda x: not any([c in x.lower() for c in ['hll','hllc','exact']]), defs_lines))
            if hll:
                extra_defs.append('#define HLL')
            if hllc:
                extra_defs.append('#define HLLC')
            if exact:
                extra_defs.append('#define EXACT')

    # Reconstruction
    pcm = int(args['pcm'])
    plm = int(args['plm'])
    ppm = int(args['ppm'])
    rsum = pcm + plm + ppm
    if rsum > 0:
        if rsum > 1:
            print('Can only have one of PCM, PLM, PPM defined!')
            exit()
        else:
            # Check for already defined Riemann solver
            defs_lines = list(filter(lambda x: not any([c in x.lower() for c in ['pcm','plm','ppm']]), defs_lines))
            if pcm:
                extra_defs.append('#define PCM')
            if plm:
                extra_defs.append('#define PLM')
            if ppm:
                extra_defs.append('#define PPM')


    add_extra_def('viscosity',args,defs_lines,extra_defs)
    add_extra_def('conduction',args,defs_lines,extra_defs)
    add_extra_def('potential',args,defs_lines,extra_defs)
    add_extra_def('ctu',args,defs_lines,extra_defs)
    add_extra_def('prof',args,defs_lines,extra_defs)
    add_extra_def('silent',args,defs_lines,extra_defs)
    add_extra_def('float',args,defs_lines,extra_defs,defname='ISFLOAT')



    # Write outputs

    shutil.copy(defsfile,'src/prob.h')
    shutil.copy(initfile,'src/prob.cu')

    lines = load_par_file(parfile)
    create_par_file(lines,'src/read_pars.c')


    struct_lines = create_struct(lines)
    with open('src/prob.h','w') as f:

        f.write('\n'.join(defs_lines + extra_defs) + '\n\n')
        f.write('\n'.join(struct_lines) + '\n')

    call(['tar','-czf','{}_src.tar.gz'.format(problem_name),'src/'])

