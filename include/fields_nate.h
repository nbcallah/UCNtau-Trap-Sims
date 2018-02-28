#ifndef UCNTFIELDS_F_H
#define UCNTFIELDS_F_H

void force_(double *x_in, double *y_in, double *z_in, double *fx, double *fy, double *fz, double *totalU, double* t, double* freq);
void fieldstrength_(double *x_in, double *y_in, double *z_in, double *totalB, double* t, double* freq);
void potential_(double *x_in, double *y_in, double *z_in, double *totalU, double* t, double* freq);

#endif /* UCNTFIELDS_F_H */
