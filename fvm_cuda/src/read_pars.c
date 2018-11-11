#include "defs.h"
#include <ctype.h>
#define PRINT_DOUBLE(NAME,VAL) printf("\t%s = %lg\n",NAME,VAL)
#define PRINT_INT(NAME,VAL) printf("\t%s = %d\n",NAME,VAL)
#define PRINT_STR(NAME,VAL) printf("\t%s = %s\n",NAME,VAL)
#define FPRINT_DOUBLE(F,NAME,VAL) fprintf(f,"%s = %lg\n",NAME,VAL)
#define FPRINT_INT(F,NAME,VAL) fprintf(f,"%s = %d\n",NAME,VAL)
#define FPRINT_STR(F,NAME,VAL) fprintf(f,"%s = %s\n",NAME,VAL)
void set_var(char *name,int int_val, double double_val, int bool_val, char *str_val, Parameters *params) {
    if (strcmp(name,"nout")==0)  { params->nout= int_val; PRINT_INT(name,int_val); }
    
    else if (strcmp(name,"nx1")==0)  { params->nx1= int_val; PRINT_INT(name,int_val); }
    
    else if (strcmp(name,"nx3")==0)  { params->nx3= int_val; PRINT_INT(name,int_val); }
    
    else if (strcmp(name,"nx2")==0)  { params->nx2= int_val; PRINT_INT(name,int_val); }
    
    else if (strcmp(name,"nscalars")==0)  { params->nscalars= int_val; PRINT_INT(name,int_val); }
        else if (strcmp(name,"cfl")==0)  { params->cfl= double_val; PRINT_DOUBLE(name,double_val); }
    
    else if (strcmp(name,"x1_max")==0)  { params->x1_max= double_val; PRINT_DOUBLE(name,double_val); }
    
    else if (strcmp(name,"x2_max")==0)  { params->x2_max= double_val; PRINT_DOUBLE(name,double_val); }
    
    else if (strcmp(name,"x2_min")==0)  { params->x2_min= double_val; PRINT_DOUBLE(name,double_val); }
    
    else if (strcmp(name,"x3_min")==0)  { params->x3_min= double_val; PRINT_DOUBLE(name,double_val); }
    
    else if (strcmp(name,"x3_max")==0)  { params->x3_max= double_val; PRINT_DOUBLE(name,double_val); }
    
    else if (strcmp(name,"x1_min")==0)  { params->x1_min= double_val; PRINT_DOUBLE(name,double_val); }
    
    else if (strcmp(name,"tend")==0)  { params->tend= double_val; PRINT_DOUBLE(name,double_val); }
    
    else if (strcmp(name,"gamma")==0)  { params->gamma= double_val; PRINT_DOUBLE(name,double_val); }
        else if (strcmp(name,"outputname")==0)  { sprintf(params->outputname,"%s",str_val); PRINT_STR(name,str_val); }
    
return;
}

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


    