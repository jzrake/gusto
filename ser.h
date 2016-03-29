#ifndef GUSTO_SERIAL
#define GUSTO_SERIAL

struct gusto_user {
    char problem_name[1024];
    char outdir[256];
    double tmax;
    double cpi;
    int N[3];
    int k2;
    double abc[3];
    double helicity;
    double pert;
    double cfl;
    double eps;
    double nu;
    int measure_cadence;
    int analyze_cadence;
    int slice_cadence;
    int num_pspec_bin;
    int max_pspec_bin;
    int calc_initial_diff;
    int normalize_initial;
};


struct gusto_status {
    int iteration;
    int checkpoint_number;
    double time_simulation;
    double time_step;
    double time_last_checkpoint;
    double kzps;
    double velocity_energy;
    double magnetic_energy;
    double velocity_monopole;
    double magnetic_monopole;
    double velocity_helicity;
    double magnetic_helicity;
    double velocity_L2;
    double magnetic_L2;
};


int gusto_user_set_from_arg(struct gusto_user *user, char *arg);
void gusto_user_report(struct gusto_user *user);
void gusto_user_set_defaults(struct gusto_user *user);
void gusto_status_set_defaults(struct gusto_status *status);


int gusto_read_write_status(struct gusto_status *status,
			    const char *chkpt_name, char mode);
int gusto_read_write_user(struct gusto_user *user,
			  const char *chkpt_name, char mode);


#endif				// GUSTO_SERIAL
