#ifndef GUSTO_SERIAL
#define GUSTO_SERIAL

struct gusto_user {
    char outdir[256];
    double tmax;
    double cpi;
    int N[2];
    int ng[2];
    double domain[4];
    char coordinates;
    int move_cells;
    char boundary_con[256];
    char initial_mesh[256];
    char initial_data[256];
    double pressure0;
    double pressure1;
    double density0;
    double density1;
    double fourvel0;
    double fourvel1;
    double sigma;
};


struct gusto_status {
    int iteration;
    int checkpoint_number;
    double time_simulation;
    double time_step;
    double time_last_checkpoint;
    double kzps;
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
