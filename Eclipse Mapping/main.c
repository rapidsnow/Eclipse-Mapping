#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define PI 3.14159265358979323846
#define DEG2RAD(DEG) ((DEG)*((PI)/(180.0)))
#include "/opt/local/include/fitsio.h"
#include "/opt/local/include/nrutil.h"
#define NR_END 1
#define FREE_ARG char*
#define NMAX 500000
#define GET_PSUM \
for (j=0;j<ndim;j++) {\
for (sum=0.0,i=0;i<mpts;i++) sum += p[i][j]; \
psum[j]=sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

/*
Versions:
 v2.1 3_11_2013 WNA
 Fixed a lot of memory issues that valgrind was reporting. It still has some leaks, but a lot
 less. There are no more explicit errors when you run valgrind on the code either.
*/

/*
 v2.03 2_19_2013 WNA
 Small bug where the lightcurve detrending flag wasn't working as expected. I didn't realize
 that argv elements are always strings, and so I was trying to use a simple Boolean comparison
 in an inappropriate place.
 
 v2.02 2_14_2013 WNA
 Lightcurve detrending is now an option at the start. You can do this by placing a 1 in the
 argument vector. This only affects outputs of files right now, not other parameters. You
 still have to manually input 0 boxes and whatnot. Binning and visibility calculations have
 been moved to the front of the code before the while loop. This should save operational
 complexity and minimize any small bugs left in those pieces of code. The normalization
 of the entire lightcurve is now down AFTER the binning so that the highest value will
 actually be one instead of .99something.
 
 v2.01 2_6_2013 WNA
 Added a test to see if the desired window size is 0. If so, do the whole data set.
 Also, fixed the binning function so that it actually bins correctly when not including
 transits (i.e. boxes set to 0). I was also defining o_cadence in the wrong place. It
 didn't matter that much for Kepler 17 because the orbital period is 1.45 days, but for
 Sarah's the orbital period is 30 days and the cadence was off by a factor of 30.
 
 v2.0 1_4_2013 WNA
 Z-values are working!
 
 v1.52 12_27_2012 WNA
 This version passes all of the necessary functions through the code properly, but
 the Z-values have not  been tested at all yet.
 
 v1.51 12_19_2012 WNA
 This version of the code has all of the necessary functions in place, but has not
 figured out how to pass them through the amoeba and the chi-squared functions properly.
 Also, more extensive commenting was added to the functions below main.
 
 v1.5 12_18_2012 WNA
 This version of the code attempts (poorly) to implement Z-Values. However, the new (and good) thing
 is that I have introduced structs to make passing things into functions easier. There are two main
 structs. s has (mostly) static inputs whose parameters are defined from the very beginning and
 do not change. It also contains arrays, but you shouldn't be modifying the arrays too much inside of
 functions with the exception of data clipping. d has the variables that change often within the
 code.
 
 v1.4 12_5_2012 WNA
 When the input file specifies 0 boxes, the code now automatically skips binning the
 transits and defining visibilities for the boxes.
 
 v1.3 12_3_2012 WNA
 Limb darkening has been included. Everything there has been checked. We are now
 normalizing box-to-box, stripe-to-stripe, and box-to-stripe. We also now are doing
 "Sliding Datasets" wherein the dataset is chopped into "windows" of a given number
 of days. Those days are then slid along the time period of the entire dataset until
 the last element of the window is greater than the last element of the dataset. I
 also reworked the binning loop so that it is a quarter of the length. The error
 binning has been recalculated.
 
 v1.2 10_3_2012 WNA
 Fixed the box values and the eclipse values. They weren't normalize correctly and were producing visibilities that weren't what they should be.

 v1.1 9_10_2012 WNA
 This code was actually updated a while ago, but it DOES have regularization.
 
 v1.0 3_12_2012 WNA
 Fitting working, could be better. No Z-Values included, no regularization.
 
 v.95 2_22_2012 WNA
 Ecliplse Path Vis Routine should be working for the most part. Made the in transit
 binning able to go down to 1 point per bin. New visibility plot outputs.
 
 v.94 2_1_2012 WNA
 WIP: Changing the eclipse vis routine fundamentally. If it fails, go back a version.
 
 v.93 1_9_2012 WNA
 WIP: Added the eclipse vis routine to the code, but have not saved in a
 separate file yet, nor have I added in the position finding, or the eclipse
 algorithm to the actual visibility code in main.
 
 v.92 10_26_2011 WNA
 Redid the math on the transit_vis subroutine. It is now correct.
 I also did a test program to ensure that it was finding proper values.
 
 v.91 10_17_2011 WNA
 Changed the indices on the arrays to where it is now 'nsb' instead of nstripes.
 This stands for number of stripes and boxes.
 Also fixed a few bugs that arose becuase of this, making the code more stable.
 Tested on a few runs, seems to be working alright
 
 v.9 10_14_2011 WNA
 Fixed longterm bug in amoeba algorithm where p was being used in place of use_p
 Took back out nboxes and replaced with only nstripes
 
 */


//====== Structs =====================================================================================================
typedef struct {
    double e_region;
    double cc;
    double a;
    double Rp;
    double lat1;
    double lat2;
    double ldc1;
    double ldc2;
    double width;
    double lam1;
    double lam2;
    double lam3;
    double orb_per;
    double o_cadence;
    double t_cadence;
    double normFactor;
    double conversion_factor;
    double days_per_tick;
    double desired_days_per_set;
    double *brights;
    double *chi;
    double *flux;
    double *orb_phase;
    double *flux_err;
    double **visibilities;
    double **bin_flux;
    double *bin_phase;
    double *bin_error;
    long npoints;
    int nboxes;
    int nstripes;
    int nsb;
    int qq;
} static_inputs;

typedef struct {
    double rph;
    double oph;
    double *EOW;
    long *nelements;
    int start;
    int stop;
    int ld_flag;
    int box_id;
    int stripe_id;
    int *c1;
    int *count;
    int *firstelem;
    int file_count;
} dynamic_inputs;

//====== Utility Functions ============================================================================================
void printerror(int status);
void free_2d(double **array, int sizeY);
void nrerror(char error_text[]);
void check_PDC_Flag(int PDC_Flag, int *err_col, int *flux_col);
void reset_array_counters(static_inputs s, dynamic_inputs d);
void reset_arrays(static_inputs s, double *y, double **simplex, long nelements);
double calculate_model_flux(static_inputs s, double*(calculate_S_vals(static_inputs)), int timestep);
double *calculate_S_vals(static_inputs s);
double **make_2d_array(int sizeY, int sizeX);
const char *get_filename_ext(const char *filename);
double getNormVal (double **bin_flux, long npoints);

//====== Conversion Functions =========================================================================================
double phase_fix(double time);

//====== Main Program Functions =======================================================================================
double stripe_vis_find(static_inputs s, dynamic_inputs d);
double box_vis_find(static_inputs s, dynamic_inputs d);
double eclipse_path_vis(static_inputs s, dynamic_inputs d, double (*box_vis)(static_inputs, dynamic_inputs));
double chi_squared(static_inputs s, dynamic_inputs d, double (*calculate_model_flux)(static_inputs, double*(*calculate_S_vals)(static_inputs), int), double*(*calculate_S_vals)(static_inputs));
void assign_visibilities(static_inputs s, dynamic_inputs d, double box_vis_find(static_inputs, dynamic_inputs), double stripe_vis_find(static_inputs, dynamic_inputs), double eclipse_path_vis(static_inputs, dynamic_inputs, double(static_inputs, dynamic_inputs)));
void bin_data(static_inputs s, dynamic_inputs d, double phase_fix(double));
void data_clipping(static_inputs s, dynamic_inputs d);
void amoeba(double **p, double y[], int ndim, double ftol, double (*funk)(static_inputs, dynamic_inputs, double(*calculate_model_flux)(static_inputs, double*(*calculate_S_vals)(static_inputs), int), double*(*calculate_S_vals)(static_inputs)), int *nfunk, static_inputs s, dynamic_inputs d, double(*calculate_model_flux)(static_inputs, double*(*calculate_S_vals)(static_inputs), int), double*(*calculate_S_vals)(static_inputs));

int main (int argc, const char * argv[]) {
	char filename[256], output_path[256], fname[256];
	long firstrow, firstelem, npoints, bjdrefi;
	int status, hdutype, flux_col, err_col, anynul, nstripes, nboxes, PDC_Flag, funk_count, start_val, end_val, points_per_day;
	double p_rot, orb_epoch, orb_per, width, nulval, o_cadence, t_cadence, scale, temp_flux, Rp, a, b, lam1, lam2, desired_days_per_set, days_per_tick, ldc1, ldc2;
	double *y, *init_bary, *init_flux, *init_error, **simplex;
    
	fitsfile *fptr;
	FILE *data;
	FILE *in;
	FILE *model;
	FILE *bright;
    FILE *hey;
    FILE *readable_input;
    FILE *chi_out;
    FILE *trans_out;
  	FILE *vis_output;
    
    static_inputs s;
    dynamic_inputs d;
	
	status	 = funk_count = 0;
	firstrow = firstelem = 1;
	nulval	 = -999;
    
    if (argv[2]) {
        strlcpy(output_path, argv[2], 256);
    } else {
        strlcpy(output_path, ".", 256);
    }
	
	in = fopen(argv[1], "r");
	if (in == NULL) {
		printf("Error: File %s not opened correctly\n", argv[1]);
		return 1;
	}
	fscanf(in, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf", filename, &p_rot, &t_cadence, &o_cadence, &orb_epoch, &width, &orb_per, &b, &Rp, &a, &lam1, &lam2, &nstripes, &nboxes, &PDC_Flag, &desired_days_per_set, &days_per_tick, &ldc1, &ldc2);
	
    //Fill in fixed values for structs
    s.nboxes = nboxes;
    s.ldc1 = ldc1;
    s.ldc2 = ldc2;
    s.nstripes = nstripes;
    s.lam1 = lam1;
    s.lam2 = lam2;
    s.nsb = s.nboxes + s.nstripes;
    s.a = a;
    s.Rp = Rp;
    s.width = 1.085*width;
    s.qq = nboxes/nstripes;
    s.o_cadence = o_cadence;
    s.t_cadence = t_cadence;
    s.orb_per = orb_per;
    if (s.nboxes) {
        s.lat1 = (b - Rp + 1)*(PI/2);
        s.lat2 = (b + Rp + 1)*(PI/2);
    } else {
        s.lat1 = s.lat2 = 0;
    }
    s.conversion_factor = s.orb_per/p_rot;
    s.e_region = ((sin(PI/2) - sin(-PI/2)) * (s.lat2 - s.lat1 - sin(s.lat2)*cos(s.lat2) + sin(s.lat1)*cos(s.lat1)))/(2 * PI);
    s.cc = (1 - s.e_region)/s.e_region;
    
    sprintf(fname, "%s/readable_inputs.in", output_path);
    readable_input = fopen(fname, "w");
    if (readable_input == NULL) {
        printf("Error: File %s not opened correctly\n", argv[1]);
		return 1;
    }
    fprintf(readable_input, "Filename:                       %s\n", filename);
    fprintf(readable_input, "Period of rotation:             %lf\n", p_rot);
    fprintf(readable_input, "In transit binning cadence:     %lf\n", t_cadence);
    fprintf(readable_input, "Out-of-transit binning cadence: %lf\n", o_cadence);
    fprintf(readable_input, "Orbital epoch:                  %lf\n", orb_epoch);
    fprintf(readable_input, "Transit width:                  %lf\n", width);
    fprintf(readable_input, "Orbital period:                 %lf\n", orb_per);
    fprintf(readable_input, "Impact parameter:               %lf\n", b);
    fprintf(readable_input, "Rp/R*:                          %lf\n", Rp);
    fprintf(readable_input, "Orbital separation:             %lf\n", a);
    fprintf(readable_input, "Regularization for stripes:     %lf\n", lam1);
    fprintf(readable_input, "Regularization for boxes:       %lf\n", lam2);
    fprintf(readable_input, "Number of stripes:              %d\n", nstripes);
    fprintf(readable_input, "Number of boxes:                %d\n", nboxes);
    fprintf(readable_input, "PDC flag:                       %d\n", PDC_Flag);
    fprintf(readable_input, "Desired Days per set:           %lf\n", desired_days_per_set);
    fprintf(readable_input, "Desired Days per tick:          %lf\n", days_per_tick);
    fprintf(readable_input, "LD Coefficient 1:               %lf\n", ldc1);
    fprintf(readable_input, "LD Coefficient 2:               %lf\n", ldc2);
    fclose(readable_input);
    
    check_PDC_Flag(PDC_Flag, &err_col, &flux_col);
	
////// Start fits files routines //////////////////////////////////////////////////////////////////
	//Open File and move to the proper header
 if (!strncmp(get_filename_ext(filename), "fits", 6)) {
        if(fits_open_file(&fptr, filename, READONLY, &status)) {
            printerror(status);
        }
        if(fits_movabs_hdu(fptr, 2, &hdutype, &status)) {
            printerror(status);
        }
        if(fits_read_key(fptr, TLONG, "NAXIS2", &npoints, NULL, &status)) {
            printerror(status);
        }
        if(fits_read_key(fptr, TLONG, "BJDREFI", &bjdrefi, NULL, &status)) {
            printerror(status);
        }
    } else {
        data = fopen(filename, "r");
        if (data == NULL) {
            printf("Error: File %s not opened correctly\n", filename);
            return 1;
        }
        fscanf(data, "%ld", &npoints);
        fscanf(data, "%ld", &bjdrefi);
    }
    
    orb_epoch = orb_epoch - bjdrefi;
    points_per_day = (int)((24*60*60)/58.85);
    
    /////////////////////////////////////////////////////////////////////////////////////////////
	if((s.orb_phase = malloc(npoints * sizeof(s.orb_phase))) == NULL) return 5;
	if((s.flux      = malloc(npoints * sizeof(s.flux)))      == NULL) return 5;
	if((s.flux_err  = malloc(npoints * sizeof(s.flux_err)))  == NULL) return 5;
	if((s.bin_phase = malloc(npoints * sizeof(s.bin_phase))) == NULL) return 5;
	s.bin_flux      = make_2d_array((int)npoints, 2);
	if((s.bin_error = malloc(npoints * sizeof(s.bin_error))) == NULL) return 5;
    if((init_bary   = malloc(npoints * sizeof(init_bary)))   == NULL) return 5;
	if((init_flux   = malloc(npoints * sizeof(init_flux)))   == NULL) return 5;
	if((init_error  = malloc(npoints * sizeof(init_error)))  == NULL) return 5;
    if((d.c1        = malloc(sizeof(d.c1)))                  == NULL) return 5;
    if((d.count     = malloc(sizeof(d.count)))               == NULL) return 5;
    if((d.firstelem = malloc(sizeof(d.firstelem)))           == NULL) return 5;
    if((d.nelements = malloc(sizeof(d.nelements)))           == NULL) return 5;
    if((d.EOW       = malloc(sizeof(d.EOW)))                 == NULL) return 5;
    simplex         = make_2d_array(s.nsb + 1, s.nsb);
    if((y           = malloc((s.nsb + 1) * sizeof(y)))       == NULL) return 5;
    if((s.chi       = malloc(sizeof(s.chi)))                 == NULL) return 5;
    /////////////////////////////////////////////////////////////////////////////////////////////
    sprintf(fname, "%s/chi.out", output_path);
    chi_out = fopen(fname, "w");
    
    *d.nelements = npoints;
    d.file_count = 0;
    s.npoints = npoints;
    s.desired_days_per_set = desired_days_per_set/s.orb_per;
    s.days_per_tick = days_per_tick/s.orb_per;
    
    if (!strncmp(get_filename_ext(filename), "fits", 6)) {
        //Read the columns into arrays
        if(fits_read_col(fptr, TDOUBLE, 1, firstrow, firstelem, npoints, &nulval, init_bary, &anynul, &status)) {
            printerror(status);		//Read Time Column
        }
        if(fits_read_col(fptr, TDOUBLE, flux_col, firstrow, firstelem, npoints, &nulval, init_flux, &anynul, &status)) {
            printerror(status);		//Read Flux Column
        }
        if(fits_read_col(fptr, TDOUBLE, err_col , firstrow, firstelem, npoints, &nulval, init_error, &anynul, &status)) {
            printerror(status);		//Read Error Column
        }
        if(fits_close_file(fptr, &status)) {
            printerror(status);
        }
    } else {
        for (int i = 0; i < npoints; i++) {
            fscanf(data, "%lf", &init_bary[i]);
            fscanf(data, "%lf", &init_flux[i]);
            fscanf(data, "%lf", &init_error[i]);
        }
        fclose(data);
    }
////// End fits files routines ////////////////////////////////////////////////////////////////////
    
    
    
///////////////////////////////////////////////////////////////////////////
    for (int i = 0; i < npoints; i++) {
        s.orb_phase[i] = init_bary[i];
        s.flux[i] = init_flux[i];
        s.flux_err[i] = init_error[i];
    }
    //Change orb_phase to units of orbital phase
    for (int k = 0; k < npoints; k++) {
        s.orb_phase[k] = (s.orb_phase[k] - orb_epoch)/s.orb_per;	//orb_phase is now in units of phase from the given orb_epoch
    }
    s.o_cadence *= (s.orb_phase[1] - s.orb_phase[0]);
    if (s.t_cadence == 1) {
        s.t_cadence = 0;
    } else {
        s.t_cadence *= (s.orb_phase[1] - s.orb_phase[0]);
    }
    
///////////////////////////////////////////////////////////////////////////
    
////////// Start Binning //////////////////////////////////////////////////////////////////////////
    //Open the binned lightcurve output file
    sprintf(fname, "%s/binned.out", output_path);
    data = fopen(fname, "w");
    
    *d.nelements = npoints;
    bin_data(s, d, phase_fix);
    
    // Normalize highest flux value to 1 //
    s.normFactor = getNormVal(s.bin_flux, *d.count);
    for (int i = 0; i <= *d.count; i++) {
        s.bin_flux[i][0] = s.bin_flux[i][0]/s.normFactor;
        s.bin_error[i] = (s.bin_error[i]/s.normFactor);
    }
    
    //Print data to file
    for (int i = 0; i <= *d.count; i++) {
        if (argc == 4) { //Lightcurve detrending flag
            fprintf(data, "%lf %lf %lf\n", s.bin_phase[i] * s.orb_per, s.bin_flux[i][0] * s.normFactor, s.bin_error[i] * s.normFactor);
        } else {
            fprintf(data, "%lf %lf %lf\n", s.bin_phase[i], s.bin_flux[i][0], s.bin_error[i]);
        }
    }
    fclose(data);
////////// End Binning ////////////////////////////////////////////////////////////////////////////

////////// Start visibilities /////////////////////////////////////////////////////////////////////
    s.visibilities = make_2d_array(*d.count + 1, s.nsb);
    d.ld_flag = 1;
    assign_visibilities(s, d, box_vis_find, stripe_vis_find, eclipse_path_vis);
    for (int i = 0; i < *d.count; i++) {
        for (int j = 0; j < s.nstripes; j++) {
            if (s.visibilities[i][j] != s.visibilities[i][j]) {
                s.visibilities[i][j] = (s.visibilities[i-1][j] + s.visibilities[i+1][j])/2;
            }
        }
    }
    
////////// End visibilities ///////////////////////////////////////////////////////////////////////
    
    
////////// Begin visibilities Outputs /////////////////////////////////////////////////////////////
    sprintf(fname, "%s/Vis_Plots/sum_bvals.out", output_path);
    vis_output = fopen(fname, "w");
    for (int j = 0; j < *d.count; j++) {
        double sum = 0;
        for (int i = 0; i < s.nboxes; i++) {
            sum += s.visibilities[j][i];
        }
        fprintf(vis_output, "%lf %lf\n", s.bin_phase[j], sum);
    }
    fclose(vis_output);
    
    sprintf(fname, "%s/Vis_Plots/sum_zvals.out", output_path);
    vis_output = fopen(fname, "w");
    for (int j = 0; j < *d.count; j++) {
        double sum = 0;
        for (int i = s.nboxes; i < s.nsb; i++) {
            sum += s.visibilities[j][i];
        }
        fprintf(vis_output, "%lf %lf\n", s.bin_phase[j], sum);
    }
    fclose(vis_output);
    
    sprintf(fname, "%s/Vis_Plots/sum_all.out", output_path);
    vis_output = fopen(fname, "w");
    for (int j = 0; j < *d.count; j++) {
        double sum = 0;
        for (int i = s.nboxes; i < s.nsb; i++) {
            sum += s.visibilities[j][i];
        }
        for (int i = 0; i < s.nboxes; i++) {
            sum -= s.visibilities[j][i];
        }
        fprintf(vis_output, "%lf %lf\n", s.bin_phase[j], sum);
    }
    fclose(vis_output);
    
    for (int i = 0; i < s.nboxes; i++) {
        sprintf(fname, "%s/Vis_Plots/b%d.out", output_path, i);
        vis_output = fopen(fname, "w");
        for (int j = 0; j < *d.count; j++) {
            fprintf(vis_output, "%lf %lf\n", s.bin_phase[j], s.visibilities[j][i]);
        }
        fclose(vis_output);
    }
    for (int i = s.nboxes; i < s.nsb; i++) {
        sprintf(fname, "%s/Vis_Plots/z%d.out", output_path, i - nboxes);
        vis_output = fopen(fname, "w");
        for (int j = 0; j < *d.count; j++) {
            fprintf(vis_output, "%lf %lf\n", s.bin_phase[j], s.visibilities[j][i]);
        }
        fclose(vis_output);
    }
////////// End visibilities Outputs ///////////////////////////////////////////////////////////////
    
    //Set up array boundaries
    *d.firstelem = 0;
    *d.EOW = s.bin_phase[0] + s.desired_days_per_set;
    for (int i = 0; i < s.npoints; i++) {
        if (s.bin_phase[i] >= *d.EOW) {
            *d.nelements = i - *d.firstelem;
            break;
        }
    }

////////////////////////////////////
//                                //
//        Start While Loop        //
//                                //
////////////////////////////////////
    
    int nn = 0;
    char lastFlag = 0;
    while (!lastFlag) {
        if (*d.firstelem + *d.nelements > *d.count) {
            lastFlag = 1;
            *d.nelements = *d.count - *d.firstelem;
        }
        nn++;
        //Print the unbinned lightcurve to an easily plottable format
        sprintf(fname, "%s/ub_%d.out", output_path, d.file_count);
        hey = fopen(fname, "w");
        for (int i = *d.firstelem; i < *d.nelements + *d.firstelem; i++) {
            if (s.flux[i] != 0) {
                fprintf(hey, "%lf %lf\n", s.orb_phase[i], s.flux[i]/s.normFactor);
            }
        }
        fclose(hey);
        
        //Assign useful variables to structs
        d.ld_flag = 1;
        
////////// Start Amoeba ///////////////////////////////////////////////////////////////////////////
        s.brights = malloc((s.nsb) * sizeof(s.brights));
        for(int i = 0; i < s.nsb; i++) {
            s.brights[i] = 1;
        }
        scale = -0.01;
        
        //Create the simplex
        for (int i = 0; i < s.nsb + 1; i++) {
            for (int j = 0; j < s.nsb; j++) {
                simplex[i][j] = s.brights[j];
                if (i == j + 1) {
                    simplex[i][j] += scale;
                }
                
            }
        }
        for (int j = 0; j < s.nsb + 1; j++) {
            memcpy(s.brights, simplex[j], s.nsb * sizeof(s.brights));
            y[j] = chi_squared(s, d, calculate_model_flux, calculate_S_vals);
        }

        amoeba(simplex, y, s.nsb, .00000000001, chi_squared, &funk_count, s, d, calculate_model_flux, calculate_S_vals);
        printf("Done with Amoeba: %d\n", nn);
        
////////// End Amoeba /////////////////////////////////////////////////////////////////////////////
        
        
////////// Start Light Curve Outputs //////////////////////////////////////////////////////////////
        if(argc == 4) {
            sprintf(fname, "%s/model_%03d.out", output_path, d.file_count);
        } else {
            sprintf(fname, "%s/model_%d.out", output_path, d.file_count);
        }
        model = fopen(fname, "w");
        for (int i = *d.firstelem; i < *d.nelements + *d.firstelem; i++) {
            temp_flux = (*calculate_model_flux)(s, calculate_S_vals, i);
            if (argc == 4) { //Lightcurve detrending flag
                fprintf(model, "%lf %lf %lf\n", s.bin_phase[i] * orb_per + orb_epoch + bjdrefi, temp_flux * s.normFactor, s.bin_error[i] * s.normFactor);
            } else {
                fprintf(model, "%lf %lf %lf\n", s.bin_phase[i], temp_flux, s.bin_error[i]);
            }
        }
        fclose(model);
        // B-Value recovery from Z-Values
        memcpy(s.brights, (*calculate_S_vals)(s), s.nsb * sizeof(s.brights));
        
        sprintf(fname, "%s/brights_%d.out", output_path, d.file_count);
        bright = fopen(fname, "w");
        for (int i = 0; i < s.nsb; i++) {
            fprintf(bright, "%lf\n", s.brights[i]);
        }
        fclose(bright);
////////// End Light Curve Outputs ////////////////////////////////////////////////////////////////
        
////////// Start Transit Range Outputs ////////////////////////////////////////////////////////////
        sprintf(fname, "%s/trans_%d.out", output_path, d.file_count);
        trans_out = fopen(fname, "w");
        for (int ii = 1 + *d.firstelem; ii < *d.nelements + *d.firstelem - 1; ii++) {
            if (s.bin_flux[ii][1] && !s.bin_flux[ii-1][1]) {
                start_val = ii;
            }
            if (s.bin_flux[ii][1] && !s.bin_flux[ii+1][1]){
                end_val = ii;
                fprintf(trans_out, "%d,%d\n", start_val + 1, end_val + 1);
            }
        }
        fclose(trans_out);
////////// End Transit Range Outputs ///////////////////////////////////////////////////////////////
        
        d.file_count++;
        funk_count = 0;
    
        reset_arrays(s, y, simplex, *d.nelements);
        reset_array_counters(s, d);
    }

    fclose(chi_out);
    free(s.orb_phase);
    free(s.flux);
    free(s.flux_err);
    free(init_bary);
    free(init_flux);
    free(init_error);
    free(y);
    free(s.bin_phase);
    free_2d(s.bin_flux, 2);
    free(s.bin_error);
    free_2d(s.visibilities, *d.count + 1);
    free_2d(simplex, s.nsb);
    
    return 0;
}

//====== Utility Functions ============================================================================================
void printerror(int status) {
    //Error reporting for cfitsio
	if(status) {
		fits_report_error(stderr,status);
		exit(status);
	}
	return;
}
void free_2d(double **array, int sizeY) {
    //Free the 2d arrays that I make with my routine below
//    for (int i = 0; i < dimension1_max; i++) {
//        free(x[i]);
//    }
//    free(x);
    for (int i = 0; i < sizeY; i++) {
        free(array[i]);
    }
	free(array);
}
void nrerror(char error_text[]) {
    //Error reporting for numerical recipes if NMAX is exceeded
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}
void check_PDC_Flag(int PDC_Flag, int *err_col, int *flux_col) {
    /* Check the PDC Flag to see which datatype to use
      If PDC_Flag = 1, use PDCSAP data otherwise SAP */
	if (PDC_Flag) {
		*flux_col = 8;
		*err_col = 9;
	} else {
		*flux_col = 4;
		*err_col = 5;
	}
}
void reset_array_counters(static_inputs s, dynamic_inputs d) {
    for (int i = *d.firstelem; i < s.npoints; i++) {
        if (s.bin_phase[i] >= s.bin_phase[*d.firstelem] + s.days_per_tick) {
            *d.firstelem = i;
            break;
        }
    }
    for (int i = *d.firstelem; i < s.npoints; i++) {
        if (s.bin_phase[i] >= *d.EOW + s.days_per_tick) {
            *d.EOW = s.bin_phase[i];
            *d.nelements = i - *d.firstelem;
            break;
        }
    }
}
void reset_arrays(static_inputs s, double *y, double **simplex, long nelements) {
    /* This was in the main program, but I didn't like the look of it
      It could be simplified, but meh. All it does is reset the array values
      so that I don't have to malloc tons of stuff every time and we don't
      use way too much memory */
    for (int i = 0; i < s.nsb + 1; i++) {
        y[i] = '\0';
        for (int j = 0; j < s.nsb; j++) {
            simplex[i][j] = '\0';
        }
    }
}
double calculate_model_flux(static_inputs s, double*(calculate_S_vals(static_inputs)), int timestep) {
    /* This calculates the model flux for one timestep.
     It should be called from a loop within either the
     chi-squared routine or the light curve reporting
     at the end of the code. */
    
    double model_flux = 0;
    double *z_bright_array;
    
    //Set one whole array to some set of "Z" vals.
    //What this array actually contains is determined by calculate_S_vals.
    //If we only want to calculate f_mod from VzZ, we should set 0 -> nboxes = 0.
    z_bright_array = (*calculate_S_vals)(s);
    
    for(int j = 0; j < s.nboxes; j++) {
        model_flux += s.visibilities[timestep][j] * z_bright_array[j];
    }
    for (int j = s.nboxes; j < s.nsb; j++) {
        double tempVis = 0;
        for (int k = 0; k < s.qq; k++) {
            tempVis += s.visibilities[timestep][(j - s.nboxes) * s.qq + k];
        }
        model_flux += (s.visibilities[timestep][j] - tempVis) * z_bright_array[j];
    }
    free(z_bright_array);
    return model_flux;
}
double *calculate_S_vals(static_inputs s) {
    /* The only instance where we need to calculate Z-values.
     This should reduce human error. */
    double tbright;
    double *brights;
    brights = malloc(s.nsb * sizeof(brights));
    
    for (int j = 0; j < s.nboxes; j++) {
        brights[j] = s.brights[j];
    }
    for (int j = s.nboxes; j < s.nsb; j++) {
        tbright = 0;
        for (int k = 0; k < s.qq; k++) {
            tbright += s.brights[(j - s.nboxes) * s.qq + k];
        }
        tbright = tbright/s.qq;
        if (s.nboxes) {
            brights[j] = (s.brights[j] - s.e_region * tbright)/(1 - s.e_region);
        } else {
            brights[j] = (s.brights[j]);
        }
    }
    return brights;
}
double **make_2d_array(int sizeY, int sizeX) {
    double **arr2d;
    
    arr2d = (double**)malloc(sizeY * sizeof(double*));
    for (int i = 0; i < sizeY; i++) {
        arr2d[i] = (double*)malloc(sizeX * sizeof(double));
    }
	return arr2d;
}
const char *get_filename_ext(const char *filename) {
    //Used to determine if the data is fits or ascii
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}
double getNormVal (double **bin_flux, long npoints) {
    /* Get a flux normalization value for the
     entire dataset, not just each individual run */
    double normFactor = bin_flux[0][0];
    for (int i = 0; i <= npoints; i++) {
        if (bin_flux[i][0] > normFactor) {
            normFactor = bin_flux[i][0];
        }
    }
    return normFactor;
}

//====== Conversion Functions =========================================================================================
double phase_fix(double time) {
    //Effectively a modulo 1 function that works with negative numbers
	if (time >= 1) {
		time -= (int)time;
	} else if (time < 0) {
		time -= (int)time;
		time += 1;
	}
	return time;
}

//====== Main Program Functions =======================================================================================
void data_clipping(static_inputs s, dynamic_inputs d) {
    /* Not operational yet */
	double mean, MAD, sum;
    mean = MAD = sum = 0;
	for (int i = d.start; i < d.stop; i++) {
		//Get a mean value for all elements in the set
		sum += s.flux[i];
		if (i == d.stop - 1) {
			mean = sum/(d.stop - d.start);
		}
	}
	
	for (int i = d.start; i < d.stop; i++) {
		//Find the absolute value of the differences and find a mean of those values
	 	sum = 0;
	 	if (s.flux[i] - mean > 0) {
	 		sum += (s.flux[i] - mean);
	 	} else {
	 		sum += (s.flux[i] - mean)*(-1);
	 	}
	 	if (i == d.stop - 1) {
	 		MAD = sum/(d.stop - d.start);
	 	}
	}
    
    //Okay we are good up through here. I understand. (11/8/12)
	for (int i = d.start; i < d.stop; i++) {
		//Set values 5 MADs away to be subtracted
	 	if (s.flux[i] < mean - (5 * MAD) || s.flux[i] > mean + (5 * MAD)) {
	 		s.bin_flux[(*d.count)][0] -= s.flux[i];
			s.bin_phase[(*d.count)] -= s.orb_phase[i];
			s.bin_error[(*d.count)] -= s.flux_err[i];
			*d.c1 -= 1;
	 	}
	}
    printf("%d %d %d\n", d.start, d.stop, *d.c1);
}
void amoeba(double **p, double y[], int ndim, double ftol, double (*funk)(static_inputs, dynamic_inputs, double(*calculate_model_flux)(static_inputs, double*(*calculate_S_vals)(static_inputs), int), double*(*calculate_S_vals)(static_inputs)), int *nfunk, static_inputs s, dynamic_inputs d, double(*calculate_model_flux)(static_inputs, double*(*calculate_S_vals)(static_inputs), int), double*(*calculate_S_vals)(static_inputs)) {
	double amotry(double **p, double y[], double psum[], int ndim, double (*funk)(static_inputs, dynamic_inputs, double(*calculate_model_flux)(static_inputs, double*(*calculate_S_vals)(static_inputs), int), double*(*calculate_S_vals)(static_inputs)), int ihi, double fac, static_inputs s, dynamic_inputs d, double(*calculate_model_flux)(static_inputs, double*(*calculate_S_vals)(static_inputs), int), double*(*calculate_S_vals)(static_inputs));
	int i, ihi, ilo, inhi, j, mpts = ndim + 1;
	double rtol, sum, swap, ysave, ytry, *psum;
	
    //use_p and use_y allow us to use the numerical recipes indexing scheme (start at 1 instead of 0)
    
    psum = malloc(ndim * sizeof(psum));
	*nfunk = 0;
	GET_PSUM
	for (;;) {
        //Testing the current max and min of each set of points
		ilo = 1;
		ihi = y[0] > y[1] ? (inhi = 1, 0) : (inhi = 0, 1);
		for (i = 0; i < mpts; i++) {
			if (y[i] <= y[ilo]) {
				ilo = i;
			}
			if (y[i] > y[ihi]) {
				inhi = ihi;
				ihi = i;
			} else if (y[i] > y[inhi] && i != ihi) {
				inhi = i;
			}
		}
        //Testing the convergence and breaking out if need be
		rtol = 2.0 * fabs(y[ihi] - y[ilo])/(fabs(y[ihi]) + fabs(y[ilo]));
		if (rtol < ftol) {
			SWAP(y[0], y[ilo])
			//SWITCHED p[ilo] TO use_p[ilo]
			for (i = 0; i < ndim; i++)  SWAP(p[0][i], p[ilo][i])
				break;
		}
        //Step limit (so that it doesn't run away)
		if (*nfunk >= NMAX) {
			nrerror("NMAX exceeded");
		}
		*nfunk += 2;
        //Modify whichever point
        
        //Reflection
		ytry = amotry(p, y, psum, ndim, funk, ihi, -1.0, s, d, calculate_model_flux, calculate_S_vals);
		if (ytry <= y[ilo]) {
            
            //Expansion
			ytry = amotry(p, y, psum, ndim, funk, ihi, 2.0, s, d, calculate_model_flux, calculate_S_vals);
		}	else if (ytry >= y[inhi]) {
			ysave = y[ihi];
            
            //One point contraction (shouldn't need 0 limiting)
			ytry = amotry(p, y, psum, ndim, funk, ihi, 0.5, s, d, calculate_model_flux, calculate_S_vals);
			if (ytry >= ysave) {
				for (i = 0; i < mpts; i++) {
					if (i != ilo) {
                        
                        //Multi Contraction (shouldn't need 0 limiting)
						for (j = 0; j < ndim; j++) {
							p[i][j] = psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
						}
                        memcpy(s.brights, psum, s.nsb * sizeof(s.brights));
						y[i] = (*funk)(s, d, calculate_model_flux, calculate_S_vals);
					}
				}
				*nfunk += ndim;
				GET_PSUM
			}
		} else {
			(*nfunk)--;
		}
	}
	free(psum);
}
void bin_data(static_inputs s, dynamic_inputs d, double phase_fix(double)) {
    /* This is only called once, but this way we can
     just trust that the data is binned. If we want
     to change something, we can do it here. */
    int transit_counter;
    int transit_status;
    int last_loop_transit_value;
    double end_of_bin;
    double half_transit;
    
    //Initialize the values... don't use calloc because it may have way too much memory for what we need.
    s.bin_flux[0][0] = 0;
    s.bin_phase[0] = 0;
    s.bin_error[0] = 0;
    
    //Start Finding Phase
    half_transit = (s.width/2)/s.orb_per;
    //End Finding Phase
    //Initialize the variables if this is the first pass through the "windows" loop
    //Reset all counting variables
    *(d.c1) = *(d.count) = transit_counter = 0;
    
    // Initialize end_of_bin and see if we are in or out-of-transit
    if(phase_fix(s.orb_phase[0]) < half_transit || 1 - phase_fix(s.orb_phase[0]) < half_transit) {
        end_of_bin = s.orb_phase[0] + s.t_cadence;
        transit_status = 1;
    } else {
        end_of_bin = s.orb_phase[0] + s.o_cadence;
        transit_status = 0;
    }
    
    /* Main loop */
    for(int i = 0; i < *d.nelements; i++) {
        if(phase_fix(s.orb_phase[i]) < half_transit || 1 - phase_fix(s.orb_phase[i]) < half_transit) {
            last_loop_transit_value = transit_status;
            transit_status = 1;
        } else {
            last_loop_transit_value = transit_status;
            transit_status = 0;
        }
        //Have we switched transits or is the bin over naturally?
        if (transit_status != last_loop_transit_value || s.orb_phase[i] > end_of_bin) {
            if (*d.c1 > 0) {
                //Are we binning the transits at all?
                if (s.nboxes == 0 && (last_loop_transit_value)) {
                    s.bin_flux[*d.count][0] = 0;
                    s.bin_phase[*d.count] = 0;
                    s.bin_error[*d.count] = 0;
                    
                    (*d.count)--;
                } else {
                    //Average the points and their errors in each bin
                    s.bin_flux[*d.count][0] = s.bin_flux[*d.count][0]/s.bin_error[*d.count];
                    s.bin_phase[*d.count] = s.bin_phase[*d.count] / *d.c1;
                    s.bin_error[*d.count] = sqrt(1/s.bin_error[*d.count]);
                    if (transit_status) {
                        //Increment the transit counter
                        s.bin_flux[*d.count][1] = transit_counter + 1;
                    } else {
                        //Set the transit counter to show that we are out-of-transit
                        s.bin_flux[*d.count][1] = 0;
                    }
                }
                //Reset the counting variables
                (*d.c1) = 0;
                (*d.count)++;
                s.bin_flux[*d.count][0] = 0;
                s.bin_phase[*d.count] = 0;
                s.bin_error[*d.count] = 0;
            }
            
            //Change value of end_of_bin
            if (transit_status != last_loop_transit_value) {
                if (transit_status) {
                    end_of_bin = s.orb_phase[i] + s.t_cadence;
                } else {
                    end_of_bin = s.orb_phase[i] + s.o_cadence;
                }
            } else {
                if (transit_status) {
                    end_of_bin += s.t_cadence;
                } else {
                    end_of_bin += s.o_cadence;
                }
            }
        }
        
        if(s.flux[i] > -990 && s.flux[i] != 0) {
            s.bin_flux[*d.count][0] += s.flux[i]/pow(s.flux_err[i],2);
            s.bin_phase[*d.count] += s.orb_phase[i];
            s.bin_error[*d.count] += 1/pow(s.flux_err[i],2);
            (*d.c1)++;
        }
    }
    /* Take care of averaging last points */
    s.bin_flux[*d.count][0] = s.bin_flux[*d.count][0]/s.bin_error[*d.count];
    s.bin_phase[*d.count] = s.bin_phase[*d.count]/ *d.c1;
    s.bin_error[*d.count] = sqrt(1/s.bin_error[*d.count]);
    
//    for (int i = 0; i <= *d.count; i++) {
//       s.bin_flux[i][0] = s.bin_flux[i][0]/s.normFactor;
//        s.bin_error[i][0] = (s.bin_error[i][0]/s.normFactor);
//    }
}
void assign_visibilities(static_inputs s, dynamic_inputs d, double box_vis_find(static_inputs, dynamic_inputs), double stripe_vis_find(static_inputs, dynamic_inputs), double eclipse_path_vis(static_inputs, dynamic_inputs, double(static_inputs, dynamic_inputs))) {
    double vis_norm;
    /* Where we define how the visibilities are set.
     Modifies the visibilities array to have a vis value
     for each timestep in each region. */
    
    for (int i = 0; i < *d.count; i++) {
        d.rph = s.bin_phase[i] * s.conversion_factor;
        d.oph = s.bin_phase[i];
        if (s.bin_flux[i][1]) {
            for (int j = 0; j < s.nboxes; j++) {
                d.box_id = j;
                s.visibilities[i][j] = box_vis_find(s, d);// - eclipse_path_vis(s, d, box_vis_find);
                s.nboxes *= 10;
                for (int k = 0; k < 10; k++) {
                    d.box_id = j * 10 + k;
                    s.visibilities[i][j] -= eclipse_path_vis(s, d, box_vis_find);
                }
                s.nboxes = s.nsb - s.nstripes;
            }
            for (int j = s.nboxes; j < s.nsb; j++) {
                d.stripe_id = j - s.nboxes;
                d.box_id = j - s.nboxes;
                s.nboxes = s.nstripes;
                
                s.visibilities[i][j] = stripe_vis_find(s, d);
                s.nboxes *= 20;
                //The reason for the 0 -> 20 loop is to give a high enough sample rate to the stripes to give a smooth limb darkening transit curve
                for (int k = 0; k < 20; k++) {
                    d.box_id = (j - s.nsb + s.nstripes) * 20 + k;
                    s.visibilities[i][j] -= eclipse_path_vis(s, d, box_vis_find);
                }
                s.nboxes = s.nsb - s.nstripes;
            }
        } else {
            for (int j = 0; j < s.nboxes; j++) {
                d.box_id = j;
                s.visibilities[i][j] = box_vis_find(s, d);
            }
            for (int j = s.nboxes; j < s.nsb; j++) {
                d.box_id = j - s.nboxes;
                d.stripe_id = j - s.nboxes;
                s.nboxes = s.nstripes;
                
                s.visibilities[i][j] = stripe_vis_find(s, d);
                
                s.nboxes = s.nsb - s.nstripes;
            }
        }
    }
    
    vis_norm = 0;
    //TODO: Doesn't work if starting in transit
    for (int i = s.nboxes; i < s.nsb; i++) {
        vis_norm += s.visibilities[0][i];
    }
    for (int i = 0; i < *d.count; i++) {
        for (int j = 0; j < s.nsb; j++) {
            s.visibilities[i][j] = s.visibilities[i][j]/vis_norm;
        }
    }
}
double stripe_vis_find(static_inputs s, dynamic_inputs d) {
	double value, start, stop, stripe_fraction_size, ld, phi1, phi2, lat1, lat2, mu_0, mu_1, mu_2;

	stripe_fraction_size = 1.0/(double)s.nstripes;
	start = phase_fix(d.rph + stripe_fraction_size * d.stripe_id); //Left edge of stripe, of Phi1/2PI
	stop = phase_fix(d.rph + stripe_fraction_size * (d.stripe_id + 1)); //Right edge of stripe, Phi2/2PI
    
    phi1 = 2 * PI * start - (PI/2);
    phi2 = 2 * PI * stop - (PI/2);
    
    lat1 = 0;
    lat2 = PI;
    
	/* Back Side */
	if (stop < 1 && stop > .5 && start < 1 && start > .5) {
        //sin(a) - sin(a) = 0, because it is on the back, ensure that the value is 0
        phi2 = 999;
        phi1 = 999;
	}
	if (stop > .5 && stop < 1 && start <= .5 && start >= 0 ) {
        //It is off the right side of the star, set the limit to PI/2
        phi2 = PI/2;
	}
	if (stop >= 0  && stop <= .5 && start < 1 && start > .5) {
        //Off the left side of the star
        phi1 = -PI/2;
	}
    
    value = (sin(phi2) - sin(phi1))/2;

    //Limb darkening
    if (phi1 != 999 && phi2 != 999) {
        mu_0 = (sin(phi2) - sin(phi1)) * (lat2 - lat1 - cos(lat2) * sin(lat2) + cos(lat1) * sin(lat1))/(2 * value * PI);
        mu_1 = (9 * (cos(lat1) - cos(lat2)) + cos(3 * lat2) - cos(3 * lat1)) * (phi2 - phi1 + sin(phi2) * cos(phi2) - sin(phi1) * cos(phi1))/(24 * value * PI);
        mu_2 = (12 * (lat2 - lat1) + 8 * (sin(2 * lat1) - sin(2 * lat2)) + sin(4 * lat2) - sin(4 * lat1)) * (9 * (sin(phi2) - sin(phi1)) + sin(3 * phi2) - sin(3 * phi1))/(384 * value * PI);
        ld = (1 - s.ldc1 - s.ldc2) * mu_0 + (s.ldc1 + 2 * s.ldc2) * mu_1 - s.ldc2 * mu_2;
    } else {
        ld = 0;
    }
    //Multiply the vis value by limb darkening
    value *= ld;
	return value;
}
double box_vis_find(static_inputs s, dynamic_inputs d) {
	double value, start, stop, box_fraction_size, phi1, phi2, mu_0, mu_1, mu_2, ld;
    
	box_fraction_size = 1.0/(double)s.nboxes;
	start = phase_fix(d.rph + box_fraction_size * d.box_id);
	stop = phase_fix(d.rph + box_fraction_size * (d.box_id + 1));
    
    phi1 = 2 * PI * start - (PI/2);
    phi2 = 2 * PI * stop - (PI/2);
    
	/* Back Side */
	if (stop < 1 && stop > .5 && start < 1 && start > .5) {
        //Set value to 0
        phi2 = 999;
        phi1 = 999;
	}
	if (stop > .5 && stop < 1 && start <= .5 && start >= 0 ) {
        //Off the right side
        phi2 = PI/2;
	}
	if (stop >= 0  && stop <= .5 && start < 1 && start > .5) {
        //Off the left side
        phi1 = -PI/2;
	}
    
    value = ((sin(phi2) - sin(phi1)) * (s.lat2 - s.lat1 - sin(s.lat2)*cos(s.lat2) + sin(s.lat1)*cos(s.lat1)))/(2 * PI);
    //Limb darkening
    if (d.ld_flag) {
        if (phi1 != 999 && phi2 != 999) {
            mu_0 = (sin(phi2) - sin(phi1)) * (s.lat2 - s.lat1 - cos(s.lat2) * sin(s.lat2) + cos(s.lat1) * sin(s.lat1))/(2 * value * PI);
            mu_1 = (9 * (cos(s.lat1) - cos(s.lat2)) + cos(3 * s.lat2) - cos(3 * s.lat1)) * (phi2 - phi1 + sin(phi2) * cos(phi2) - sin(phi1) * cos(phi1))/(24 * value * PI);
            mu_2 = (12 * (s.lat2 - s.lat1) + 8 * (sin(2 * s.lat1) - sin(2 * s.lat2)) + sin(4 * s.lat2) - sin(4 * s.lat1)) * (9 * (sin(phi2) - sin(phi1)) + sin(3 * phi2) - sin(3 * phi1))/(384 * value * PI);
            ld = (1 - s.ldc1 - s.ldc2) * mu_0 + (s.ldc1 + 2 * s.ldc2) * mu_1 - s.ldc2 * mu_2;
        } else {
            ld = 0;
        }
        value *= ld;
    }
	return value;
}
double eclipse_path_vis(static_inputs s, dynamic_inputs d, double (*box_vis)(static_inputs, dynamic_inputs)) {
	double rpart, value, inc;
	int c_num;
    
	inc = 0;
	//Pre-test
	double start, stop, box_fraction_size, ratio;
	box_fraction_size = 1.0/(double)s.nboxes;
	start = phase_fix(d.rph + box_fraction_size *  d.box_id);
	stop = phase_fix(d.rph + box_fraction_size * (d.box_id + 1));
	d.oph = phase_fix(d.oph);
	
    //This is how we include limb darkening for the transits
    //It is kind of "cheating"
    //Take the ration of a limb-darkened box/non-limb-darkened box
    dynamic_inputs d_no_ld;
    memcpy(&d_no_ld, &d, sizeof(d));
    d_no_ld.ld_flag = 0;
    ratio = box_vis(s, d)/box_vis(s, d_no_ld);
    
	/* Back Side */
	if (stop < 1 && stop > .5 && start < 1 && start > .5) {
		return 0;
	}
	
	//Variables
	double ly, ry, y0, test_value_left, test_value_right, theta, right_bound_test, left_bound_test;
	
	theta = (s.lat1 + s.lat2)/2;
    
	left_bound_test = (2 * PI * (d.rph +  d.box_id/(double)s.nboxes));
	right_bound_test = (2 * PI * (d.rph + (d.box_id + 1)/(double)s.nboxes));
	//Star goes from 0 to PI in these values
	
	
	if (left_bound_test > 0) {
		left_bound_test = fmod(left_bound_test, 2 * PI);
	} else {
		left_bound_test = 2 * PI - fmod(-left_bound_test, 2 * PI);
	}
    
	if (right_bound_test > 0) {
		right_bound_test = fmod(right_bound_test, 2 * PI);
	} else {
		right_bound_test = 2 * PI - fmod(-right_bound_test, 2 * PI);
	}
	//The Domain of left and right bound tests are 0 -> 2PI. The front side is 0 -> PI. The backside is PI -> 2PI
	
	
	ly = sin(2 * PI * (d.rph +  d.box_id/(double)s.nboxes) - PI/2) * sin(theta);
    
	ry = sin(2 * PI * (d.rph + (d.box_id + 1)/(double)s.nboxes) - PI/2) * sin(theta);
	//Front side of ry and ly is -PI/2 -> PI/2. Backside is -PI -> -PI/2 and PI/2 -> PI
	
	y0 = s.a * sin(d.oph * 2 * PI);
	
    
	if (left_bound_test > PI && left_bound_test < (2 * PI) && right_bound_test >= 0 && right_bound_test <= PI) {
		ly = sin(-PI/2) * sin(theta);
	}
	if (right_bound_test > PI && right_bound_test < (2 * PI) && left_bound_test >= 0 && left_bound_test <= PI) {
		ry = sin(PI/2) * sin(theta);
	}
    
	test_value_left = y0 - ly; //Represents distance from left stripe to center of planet
	test_value_right = y0 - ry; //Represents distance from right stripe to center of planet
	//For both of the above values, positive is when the stripe side is on the left of the planet:
	//   |--->. The ---> is positive when "|" is a side of a stripe and "." is the center of the planet.
	
	//Case 1: Planet is completely out of box to the right
	if(test_value_left > s.Rp && test_value_right > s.Rp) {
		value = 0; //Nothing changes
		c_num = 1;
	}
	//Case 2: Planet is completely out of box to the left
	if (test_value_left < -s.Rp && test_value_right < -s.Rp) {
		value = 0; //Nothing changes
		c_num = 2;
	}
	
	//===================================================//
	//Case 3: Planet is completely contained in a box
	if(test_value_left >= s.Rp && test_value_right <= -s.Rp) {
		value = PI * pow(s.Rp, 2);
		c_num = 3;
	}
	//Case 4: Planet is partially off the right side of the box
	if(test_value_left >= s.Rp && test_value_right > -s.Rp && test_value_right <= s.Rp) {
		rpart = fabs(test_value_right);
		//Case 4a: Center of planet is inside box
		if (test_value_right <= 0) {
			value = PI * pow(s.Rp, 2) - (pow(s.Rp, 2)*(acos(rpart/s.Rp)) - rpart*sqrt(pow(s.Rp, 2) - pow(rpart, 2)));
			c_num = 41;
		}
		//Case 4b: Center of planet is outside of box to the right
		if (test_value_right > 0) {
			value = pow(s.Rp, 2)*(acos(rpart/s.Rp)) - rpart*sqrt(pow(s.Rp, 2) - pow(rpart, 2));
			c_num = 42;
		}
	}
	//Case 5: Planet is partially off the left side of the box
	if(test_value_left >= -s.Rp && test_value_left < s.Rp && test_value_right <= -s.Rp) {
		rpart = fabs(test_value_left);
		//Case 5a: Center of planet is inside box
		if (test_value_left >= 0) {
			value = PI * pow(s.Rp, 2) - (pow(s.Rp, 2)*(acos(rpart/s.Rp)) - rpart*sqrt(pow(s.Rp, 2) - pow(rpart, 2)));
			c_num = 51;
		}
		//Case 5b: Center of planet is outside of box to the left
		if (test_value_left < 0) {
			value = pow(s.Rp, 2)*(acos(rpart/s.Rp)) - rpart*sqrt(pow(s.Rp, 2) - pow(rpart, 2));
			c_num = 52;
		}
	}
	//Case 6: Planet is contained partially in a box, but neither limb is in the box
	if(test_value_left < s.Rp && test_value_left > -s.Rp && test_value_right > -s.Rp && test_value_right < s.Rp) {
		//Case 6a: Planet center is in box
		if (test_value_left >= 0 && test_value_right <= 0) {
			rpart = fabs(test_value_left);
			value = PI * pow(s.Rp, 2) - (pow(s.Rp, 2)*(acos(rpart/s.Rp)) - rpart*sqrt(pow(s.Rp, 2) - pow(rpart, 2)));
			rpart = fabs(test_value_right);
			value = value - (pow(s.Rp, 2)*(acos(rpart/s.Rp)) - rpart*sqrt(pow(s.Rp, 2) - pow(rpart, 2)));
			c_num = 61;
		}
		//Case 6b: Planet center is to the left of the box
		if (test_value_left < 0) {
			rpart = fabs(test_value_left);
			value = pow(s.Rp, 2)*(acos(rpart/s.Rp)) - rpart*sqrt(pow(s.Rp, 2) - pow(rpart, 2));
			rpart = fabs(test_value_right);
			value = value - (pow(s.Rp, 2)*(acos(rpart/s.Rp)) - rpart*sqrt(pow(s.Rp, 2) - pow(rpart, 2)));
			c_num = 62;
		}
		//Case 6c: Planet center is to the right of the box
		if (test_value_right > 0) {
			rpart = fabs(test_value_right);
			value = pow(s.Rp, 2)*(acos(rpart/s.Rp)) - rpart*sqrt(pow(s.Rp, 2) - pow(rpart, 2));
			rpart = fabs(test_value_left);
			value = value - (pow(s.Rp, 2)*(acos(rpart/s.Rp)) - rpart*sqrt(pow(s.Rp, 2) - pow(rpart, 2)));
			c_num = 63;
		}
	}
    //Normalize the result by dividing by PI
    value = value/PI;
    //Multiply by the average limb darkening ratio in this region at this given timestep
    value *= ratio;
	return value;
}
double chi_squared(static_inputs s, dynamic_inputs d, double (*calculate_model_flux)(static_inputs, double*(*calculate_S_vals)(static_inputs), int), double*(*calculate_S_vals)(static_inputs)) {
    /* The function that we are working to minimize in the amoeba
     algorithm. Right now it looks really complicated because of
     the chi squared reporting for regularizatin. In reality, we
     don't need those if statements and it could be very simple.
     
     Returns: chi-squared value for the given brightness guesses */
    
	double temp_flux, weight, tflux1, tflux2;
    
    *s.chi = 0;
	for (int i = *d.firstelem; i < *d.firstelem + *d.nelements; i++) {
		temp_flux = (*calculate_model_flux)(s, calculate_S_vals, i);
        
		//Weight the chi squared so that the transits are more important
		if (s.bin_flux[i][1]) {
			weight = 1;
            *s.chi += pow((s.bin_flux[i][0] - temp_flux),2)/pow(s.bin_error[i], 2) * weight;
            if (s.lam1 + s.lam2) {
                //Regularize the light curve, first for stripes
                for (int j = s.nboxes; j < s.nsb; j++) {
                    if (j == s.nsb - 1) {
                        tflux1 = s.brights[j];
                        tflux2 = s.brights[s.nboxes];
                    } else {
                        tflux1 = s.brights[j];
                        tflux2 = s.brights[j+1];
                    }
                    *s.chi += s.lam1 * pow(tflux2 - tflux1, 2) * weight;
                }
                
                 //Then for boxes
                 for (int j = 0; j < s.nboxes; j++) {
                     if (j == s.nboxes - 1) {
                         tflux1 = s.brights[j];
                         tflux2 = s.brights[0];
                     } else {
                         tflux1 = s.brights[j];
                         tflux2 = s.brights[j+1];
                     }
                     *s.chi += s.lam2 * pow(tflux2 - tflux1, 2) * weight;
                 }
                
                //And for boxes next to stripes
                for (int j = s.nboxes; j < s.nsb; j++) {
                    for (int k = 0; k < s.qq; k++) {
                        tflux1 = s.brights[j];
                        tflux2 = s.brights[(j - s.nboxes) * s.qq + k];
                        *s.chi += s.lam2 * pow(tflux2 - tflux1, 2) * weight;
                    }
                }
            }
		} else {
			weight = 1;
            *s.chi += pow((s.bin_flux[i][0] - temp_flux),2)/pow(s.bin_error[i], 2) * weight;
            if (s.lam1 + s.lam2) {
                //Regularize the light curve, first for stripes
                for (int j = s.nboxes; j < s.nsb; j++) {
                    if (j == s.nsb - 1) {
                        tflux1 = s.brights[j];
                        tflux2 = s.brights[s.nboxes];
                    } else {
                        tflux1 = s.brights[j];
                        tflux2 = s.brights[j+1];
                    }
                    *s.chi += s.lam1 * pow(tflux2 - tflux1, 2) * weight;
                }
                
                 //Then for boxes
                 for (int j = 0; j < s.nboxes; j++) {
                     if (j == s.nboxes - 1) {
                         tflux1 = s.brights[j];
                         tflux2 = s.brights[0];
                     } else {
                         tflux1 = s.brights[j];
                         tflux2 = s.brights[j+1];
                     }
                     *s.chi += s.lam2 * pow(tflux2 - tflux1, 2) * weight;
                 }
                
                //And for boxes next to stripes
                for (int j = s.nboxes; j < s.nsb; j++) {
                    for (int k = 0; k < s.qq; k++) {
                        tflux1 = s.brights[j];
                        tflux2 = s.brights[(j - s.nboxes) * s.qq + k];
                        *s.chi += s.lam2 * pow(tflux2 - tflux1, 2) * weight;
                    }
                }
            }
		}
	}
	return *s.chi;
}
double amotry(double **p, double y[], double psum[], int ndim, double (*funk)(static_inputs, dynamic_inputs, double(*calculate_model_flux)(static_inputs, double*(*calculate_S_vals)(static_inputs), int), double*(*calculate_S_vals)(static_inputs)), int ihi, double fac, static_inputs s, dynamic_inputs d, double(*calculate_model_flux)(static_inputs, double*(*calculate_S_vals)(static_inputs), int), double*(*calculate_S_vals)(static_inputs)) {
	int j;
	double fac1, fac2, ytry, *ptry;
	
	ptry = malloc(ndim * sizeof(ptry));
	fac1 = (1.0 - fac)/ndim;
	fac2 = fac1 - fac;
	for (j = 0; j < ndim; j++) {
		ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;
	}
    memcpy(s.brights, ptry, s.nsb * sizeof(s.brights));
	ytry = (*funk)(s, d, calculate_model_flux, calculate_S_vals);
	if (ytry < y[ihi]) {
		y[ihi] = ytry;
		for (j = 0; j < ndim; j++) {
			psum[j] += ptry[j] - p[ihi][j];
			p[ihi][j] = ptry[j];
		}
	}
	free(ptry);
	return ytry;
}

