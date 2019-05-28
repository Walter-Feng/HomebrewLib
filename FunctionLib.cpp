#include "FunctionLib.h"

///////////////////////// Definition of the structures, mainly used as parameters for functions ///////////////////////


///////////Meaning of Words:////////////////////
//
//  coef ---> coefficient
//  transl ---> translation
//  stddev ---> standard deviation
//  magnif ---> magnification of the variable, 
//  e.g.: sin(2x) where 2 is the magnification
//
//  in all functions the combination of magnification 
//  and translation is interpreted as
//  f( magnif( x - transl) ).
//
////////////////////////////////////////////////



// Function for polynomials
double PolyFunction(double x, void * param)
{
    double sum;
    int i;

    PolyFunc * tempparam = (PolyFunc *) param;
    sum = 0.0;

    for(i=0;i<tempparam->length;i++)
        sum += (tempparam->coef)[i] * gsl_pow_int((tempparam->magnif)*(x - tempparam->transl), i);

    return sum;
}

// Function for Gaussian Function
double GaussianFunction(double x,void * param)
{
    GaussianFunc * tempparam = (GaussianFunc *) param;

    double result;
    
    result = tempparam->coef * exp(-gsl_pow_2((x-tempparam->mean)/tempparam->stddev)/2);

    return result;
}

//Function for sinusoidal series
double SinCosFunction(double x,void * param)
{
    SinCosFunc * tempparam = (SinCosFunc *) param;

    double sum;
    
    int i;

    sum = 0.0;
    
    for(i=0;i<tempparam->length;i++){
        sum += (tempparam->sincoef)[i] * sin((double) i *(tempparam->magnif * (x - tempparam->transl)));

        sum += (tempparam->coscoef)[i] * cos((double) i *(tempparam->magnif * (x - tempparam->transl)));
    }

    return sum;
}

double Exponential(double x,void * param)
{
    Exp * tempparam = (Exp *) param;

    double result;

    result = tempparam->coef * exp((x-tempparam->transl)*tempparam->magnif);

    return result;
}

double Logarithm(double x,void * param)
{
    Log * tempparam = (Log *) param;
    double result;

    result = tempparam->coef * log((x-tempparam->coef)/tempparam->magnif) / log(tempparam->base);
    
    return result;
}

double Rectangular(double x, void * param)
{
    Rect * tempparam = (Rect *) param;

    double result;
    if(tempparam->a <= x && x<= tempparam->b)
        return tempparam->height;
    else return 0;
}
double Combination(int n,int r)
{
    return gsl_sf_fact(n)/gsl_sf_fact(r)/gsl_sf_fact(n-r);
}

void LegendrePCoefficient(PolyFunc * target, LegendreP * param)
{
    int i;
    double temp;

    i = 0;

    target->length = param->l + 1;

    if(target->coef != NULL)
        free (target->coef);

    target->coef = new double[param->l + 1];

    for(i=0;i<param->l+1;i++)
        (target->coef[i]) = 0;

    i=0;
   while(i<(param->l)/2+1){
        temp = (param->coef) * gsl_pow_int(-1.0,i) * gsl_sf_fact(2*param->l - 2 * i)/gsl_sf_fact(i)/gsl_sf_fact(param->l - i)/gsl_sf_fact(param->l-2*i)/gsl_pow_int(2.0,param->l);

        (target->coef)[param->l - 2*i] = temp;

        i++;
    }

    target->transl = param->transl;
    target->magnif = param->magnif;
}

void gsl_vector_poly_convert(PolyFunc * target,gsl_vector * source, int vectorlength)
{
    int i;

    if(target->coef != NULL)
        delete target->coef;

    target->coef = new double[vectorlength];

    for(i=0;i<vectorlength;i++)
        (target->coef)[i] = gsl_vector_get(source,i);

    target->length = vectorlength;
    target->transl = 0;
    target->magnif = 1;
}

void PolynomialTranslationCoef(gsl_vector * target,PolyFunc * param,int vectorlength)
{
    gsl_vector * tempvector;
    
    int i,j;

    double temp1,temp2;
    double mag;

    gsl_vector_set_zero(target);

    tempvector = gsl_vector_calloc(vectorlength);

    for(i=0;i<param->length;i++)
        gsl_vector_set(tempvector,i,(param->coef)[i]);

    for(i=0;i<param->length;i++){
        mag = gsl_pow_int(param->magnif,i);
        for(j=0;j<=i;j++){
            temp1 = gsl_vector_get(tempvector,i) * mag * Combination(i,j) * gsl_pow_int(-param->transl,i-j);

            temp2 = gsl_vector_get(target,j);

            gsl_vector_set(target,j,temp1+temp2);
        }
    }
    gsl_vector_free(tempvector);
}



void OrthoPolyCoefGen(gsl_matrix * target,int rows,int columns,double a,double b)
{
    gsl_vector * vectortemp;

    int i;

    LegendreP generator;
    PolyFunc converted;

    converted.coef = NULL;


    generator.magnif = 2.0/(b-a);
    generator.transl = (b+a)/2.0;
    generator.coef = 1.0;

    vectortemp = gsl_vector_calloc(rows);

    gsl_matrix_set_zero(target);



    for(i=0;i<columns;i++){
        generator.l = i;
        generator.coef = sqrt(2.0/(b-a)) * sqrt((2.0*(double) i + 1.0)/2.0);
        LegendrePCoefficient(&converted,&generator);
        PolynomialTranslationCoef(vectortemp,&converted,rows);
        gsl_matrix_set_col(target,i,vectortemp);
    }

    gsl_vector_free(vectortemp);
}

double IDMapping(double x, int ID, void * param)
{
    double result;

    switch(ID)
    {
        case 1: result = PolyFunction(x,param);break;
        case 2: result = GaussianFunction(x,param);break;
        case 3: result = SinCosFunction(x,param);break;
        case 4: result = Exponential(x,param);
        case 5: result = Logarithm(x,param);break;
        case 7: result = Rectangular(x,param);break;
    }

    return result;
}

double FunctionMultiply(double x, void * HeadMulMap)
{
    MulMap * tempparam;
    double result;

    if(HeadMulMap == NULL)
        return 0.0;

    result = 1.0;

    tempparam = (MulMap *) HeadMulMap;

    if(tempparam->flag == 1)
        result = result * IDMapping(x,tempparam->ID,tempparam->target);
    else result = result / IDMapping(x,tempparam->ID,tempparam->target);

    while(tempparam->NEXT != NULL){
        tempparam = (MulMap *) tempparam->NEXT;
        if(tempparam->flag == 1)
            result = result * IDMapping(x,tempparam->ID,tempparam->target);
        else result = result / IDMapping(x,tempparam->ID,tempparam->target);
    }
    
    return result;
}

double FunctionSum(double x,void * summap)
{
    SumMap * summaptemp;
    MulMap * mulmaptemp;
    double result;

    summaptemp = (SumMap *) summap;
    
    mulmaptemp = (MulMap *) summaptemp->MulMaphead;

    result = 0;

    if(summaptemp->flag == 1)
            result = result + FunctionMultiply(x,mulmaptemp);
        else result = result - FunctionMultiply(x,mulmaptemp);

    while(summaptemp->NEXT != NULL){
        summaptemp = (SumMap *) summaptemp->NEXT;
        mulmaptemp = (MulMap *) summaptemp->MulMaphead;
        if(summaptemp->flag == 1)
            result = result + FunctionMultiply(x,mulmaptemp);
        else result = result - FunctionMultiply(x,mulmaptemp);

    }

    return result;
}

void FuncStructDelete(void* target, int ID)
{
    switch(ID)
    {
        case 1: delete ((PolyFunc *) target)->coef; delete (PolyFunc *) target; break;
        case 2: delete (GaussianFunc *) target; break;
        case 3: delete ((SinCosFunc *) target)->sincoef; delete ((SinCosFunc *) target)->coscoef; delete (SinCosFunc *) target; break;
        case 4: delete (Exp *) target; break;
        case 5: delete (Log *) target; break;
        case 7: delete (Rect *) target; break;
    } 
}

void MulMapDelete(void * target)
{
    MulMap * mulmapdel;
    MulMap * mulmapnext;

    mulmapdel = (MulMap *) target;
    

    while(mulmapdel->NEXT != NULL){
        mulmapnext = (MulMap *) mulmapdel->NEXT;
        FuncStructDelete(mulmapdel->target,mulmapdel->ID);
        delete mulmapdel;
        mulmapdel = mulmapnext;
    }

    FuncStructDelete(mulmapdel->target,mulmapdel->ID);

    delete mulmapdel;

    ((MulMap *) target)->NEXT = NULL;
}

void SumMapDelete(void * target)
{
    SumMap * summapdel;
    SumMap * summapnext;

    summapdel = (SumMap *) target;

    while(summapdel->NEXT != NULL){
        summapnext = (SumMap*) summapdel->NEXT;
        MulMapDelete(summapdel->MulMaphead);
        delete summapdel;
        summapdel = summapnext;
    }

    MulMapDelete(summapdel->MulMaphead);

    delete summapdel;

    ((SumMap *) target)->NEXT = NULL;
}

void FuncStructConstruct(FILE * inputfile,void ** target, int ID)
{
    int i,length;
    double translation;
    double magnification;
    double * coefficient;
    
    coefficient = NULL;
    
    switch(ID)
    {
        case 1:
            (* target) =(void *) (new PolyFunc);
            fscanf(inputfile,"%d",&length);
            coefficient = new double[length];
            for(i=0;i<length;i++){
                fscanf(inputfile,"%lf",&(coefficient[i]));
            }
            fscanf(inputfile,"%lf",&translation);
            fscanf(inputfile,"%lf",&magnification);
            ((PolyFunc *) (* target))->length = length;
            ((PolyFunc *) (* target))->coef = coefficient;
            ((PolyFunc *) (* target))->transl = translation;
            ((PolyFunc *) (* target))->magnif = magnification;
            break;
            
        case 2:
            (* target) = (void *) (new GaussianFunc);
            fscanf(inputfile,"%lf",&((GaussianFunc *)(* target))->mean);
            fscanf(inputfile,"%lf",&((GaussianFunc *)(* target))->stddev);
            fscanf(inputfile,"%lf",&((GaussianFunc *)(* target))->coef);
            break;
            
        case 3:
            (* target) = (void *) (new SinCosFunc);
            fscanf(inputfile,"%d",&length);
            coefficient = new double[length];
            for(i=0;i<length;i++)
                fscanf(inputfile,"%lf",&(coefficient[i]));
            ((SinCosFunc *) (* target))->sincoef = coefficient;
            coefficient = new double[length];
            for(i=0;i<length;i++)
                fscanf(inputfile,"%lf",&(coefficient[i]));
            ((SinCosFunc *) (* target))->coscoef = coefficient;
            fscanf(inputfile,"%lf",&((SinCosFunc *)(* target))->transl);
            fscanf(inputfile,"%lf",&((SinCosFunc *)(* target))->magnif);
            break;
            
        case 4:
            (* target) = (void *) (new Exp);
            fscanf(inputfile,"%lf",&((Exp *)(* target))->coef);
            fscanf(inputfile,"%lf",&((Exp *)(* target))->transl);
            fscanf(inputfile,"%lf",&((Exp *)(* target))->magnif);
            break;
            
        case 5:
            (* target) = (void *) (new Log);
            fscanf(inputfile,"%lf",&((Log *)(* target))->coef);
            fscanf(inputfile,"%lf",&((Log *)(* target))->transl);
            fscanf(inputfile,"%lf",&((Log *)(* target))->base);
            fscanf(inputfile,"%lf",&((Log *)(* target))->magnif);
            break;
            
        case 6:
            (* target) = (void *) (new LegendreP);
            fscanf(inputfile,"%d",&((LegendreP *)(* target))->l);
            fscanf(inputfile,"%lf",&((LegendreP *)(* target))->transl);
            fscanf(inputfile,"%lf",&((LegendreP *)(* target))->magnif);
            break;

        case 7:
            (* target) = (void *) (new Rect);
            fscanf(inputfile,"%lf",&((Rect *)(* target))->a);
            fscanf(inputfile,"%lf",&((Rect *)(* target))->b);
            fscanf(inputfile,"%lf",&((Rect *)(* target))->height);
    }
}

void MulMapConstruct(FILE * inputfile, MulMap * headpointer)
{
    char tempstr[20];
    int ID;
    int flag;
    void * funcstruct = NULL;
    MulMap * mulmaptemp1;
    MulMap * end;
    
    mulmaptemp1 = headpointer;
    
    while(fscanf(inputfile,"%s",tempstr)){
        if(strcmp(tempstr,"BEGINMULMAP")==0)
		continue;
	else if(strcmp(tempstr,"ENDMULMAP")==0)
            break;
        else if(strcmp(tempstr,"FUNCCALL")==0){
            fscanf(inputfile,"%s",tempstr);
            fscanf(inputfile,"%d",&flag);
            fscanf(inputfile,"%s",tempstr);
            fscanf(inputfile,"%d",&ID);
            FuncStructConstruct(inputfile,&funcstruct,ID);
            
            mulmaptemp1->target = (void *) funcstruct;
            mulmaptemp1->flag = flag;
            mulmaptemp1->ID = ID;
            mulmaptemp1->NEXT = (void *) new MulMap;
            end = mulmaptemp1;
            mulmaptemp1 = (MulMap *) mulmaptemp1->NEXT;
        }
    }
    delete (MulMap *) end->NEXT;
    end->NEXT = NULL;
    
}

void SumMapConstruct(FILE * inputfile, SumMap * headpointer)
{
    MulMap * mulmaptemp = new MulMap;
    SumMap * summaptemp1;
    SumMap * summapend;
    
    int flag;
    
    char tempstr[20];
    
    summaptemp1 = headpointer;
    
    while(fscanf(inputfile,"%s",tempstr)){
        if(strcmp(tempstr,"BEGINSUMMAP")==0){
            fscanf(inputfile,"%s",tempstr);
            fscanf(inputfile,"%d",&flag);
            summaptemp1->flag = flag;
        }
        
        else if(strcmp(tempstr,"MULMAPCALL")==0)
        {
            MulMapConstruct(inputfile,mulmaptemp);
            summaptemp1->MulMaphead = (void *) mulmaptemp;
            summaptemp1->NEXT = (void *) new SumMap;
            summapend = summaptemp1;
            summaptemp1 = (SumMap *) summaptemp1->NEXT;
            mulmaptemp = new MulMap;
        }
        else if(strcmp(tempstr,"ENDSUMMAP")==0)
            break;
        
    }
    delete (SumMap *) summapend->NEXT;
    summapend->NEXT = NULL;
    delete mulmaptemp;
    
}

