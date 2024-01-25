#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<stdint.h>
#include<unistd.h>
#include<sys/wait.h>
#include<mpfr.h>
#include<pthread.h>

struct ctx {
  mpq_t L, Lc;
  mpq_t X, Xc;
  mpq_t K, Kc;
  mpq_t M, Mc;
  mpq_t T1;
  mpq_t T2;
  mpq_t T3;
  int i;

  mpfr_t pi;
};

void *adv_L(void *p) {
  struct ctx *c = p;
  mpq_add(c->L, c->L, c->Lc);
  mpq_canonicalize(c->L);
  return NULL;
}

void *adv_X(void *p) {
  struct ctx *c = p;
  mpq_mul(c->X, c->X, c->Xc);
  mpq_canonicalize(c->X);
  return NULL;
}

void *adv_M(void *p) {
  struct ctx *c = p;
  mpq_mul(c->T1, c->K, c->K);
  mpq_mul(c->T1, c->T1, c->K);
  mpq_mul(c->T2, c->Mc, c->K);
  mpq_sub(c->T1, c->T1, c->T2);
  mpq_set_ui(c->T2, c->i + 1, 1);
  mpq_div(c->T1, c->T1, c->T2);
  mpq_div(c->T1, c->T1, c->T2);
  mpq_div(c->T1, c->T1, c->T2);
  mpq_mul(c->M, c->M, c->T1);
  mpq_canonicalize(c->M);
  return NULL;
}

void *adv_K(void *p) {
  struct ctx *c = p;
  mpq_add(c->K, c->K, c->Kc);
  mpq_canonicalize(c->K);
  return NULL;
}

void *calc_sum_elem(void *p) {
  struct ctx *c = p;
  mpq_mul(c->T3, c->M, c->L);
  mpq_div(c->T3, c->T3, c->X);
  mpq_canonicalize(c->T3);
  return NULL;
}

void *add_to_sum(void *p) {
  struct ctx *c = p;
  mpfr_add_q(c->pi, c->pi, c->T3, MPFR_RNDN);
  return NULL;
}

void calc_pi(mpfr_t pi, int precision) {
  pthread_t thL, thX, thM, thK, thS, thPi;
  int bits = (precision * 11) / 3;

  MPFR_DECL_INIT(C, bits);
  mpfr_sqrt_ui(C, 10005, MPFR_RNDN);
  mpfr_mul_ui(C, C, 426880, MPFR_RNDN);

  struct ctx c;

  mpfr_init2(c.pi, bits);
  mpfr_set_ui(c.pi, 0, MPFR_RNDN);

  mpq_init(c.L);
  mpq_init(c.X);
  mpq_init(c.K);
  mpq_init(c.M);
  mpq_init(c.Lc);
  mpq_init(c.Xc);
  mpq_init(c.Kc);
  mpq_init(c.Mc);
  mpq_init(c.T1);
  mpq_init(c.T2);
  mpq_init(c.T3);

  mpq_set_ui(c.Lc, 545140134, 1);
  
  mpq_set_ui(c.Xc, 640320, 1);
  mpq_mul(c.T1, c.Xc, c.Xc);
  mpq_mul(c.Xc, c.Xc, c.T1);
  mpq_neg(c.Xc, c.Xc);

  mpq_set_ui(c.Kc, 12, 1);
  mpq_set_ui(c.Mc, 16, 1);

  mpq_set_ui(c.L, 13591409, 1);
  mpq_set_ui(c.X, 1, 1);
  mpq_set_si(c.K, -6, 1);
  mpq_set_ui(c.M, 1, 1);

  for(int i = 0; i < precision / 14; i++) {
    c.i = i;

    pthread_create(&thS, NULL, calc_sum_elem, &c);
    pthread_create(&thK, NULL, adv_K, &c);
#ifdef DEBUG
    fprintf(stderr, "Iter#%d of %d size: L:%ld, X:%ld, K:%ld, M:%ld, T1:%ld, T2:%ld, T3:%ld\n",
        i, precision / 14,
        mpz_size(mpq_numref(c.L)) + mpz_size(mpq_denref(c.L)),
        mpz_size(mpq_numref(c.X)) + mpz_size(mpq_denref(c.X)),
        mpz_size(mpq_numref(c.K)) + mpz_size(mpq_denref(c.K)),
        mpz_size(mpq_numref(c.M)) + mpz_size(mpq_denref(c.M)),
        mpz_size(mpq_numref(c.T1)) + mpz_size(mpq_denref(c.T1)),
        mpz_size(mpq_numref(c.T2)) + mpz_size(mpq_denref(c.T2)),
        mpz_size(mpq_numref(c.T3)) + mpz_size(mpq_denref(c.T3)));
#endif /* DEBUG */
    pthread_join(thS, NULL);
    pthread_create(&thPi, NULL, add_to_sum, &c);

    pthread_create(&thL, NULL, adv_L, &c);
    pthread_create(&thX, NULL, adv_X, &c);
    pthread_join(thK, NULL);
    pthread_create(&thM, NULL, adv_M, &c);

    pthread_join(thL, NULL);
    pthread_join(thX, NULL);
    pthread_join(thM, NULL);
    pthread_join(thPi, NULL);
  }

  calc_sum_elem(&c);
  add_to_sum(&c);

  mpq_clear(c.L);
  mpq_clear(c.X);
  mpq_clear(c.K);
  mpq_clear(c.M);
  mpq_clear(c.Lc);
  mpq_clear(c.Xc);
  mpq_clear(c.Kc);
  mpq_clear(c.Mc);
  mpq_clear(c.T1);
  mpq_clear(c.T2);
  mpq_clear(c.T3);

  mpfr_ui_div(c.pi, 1, c.pi, MPFR_RNDN);
  mpfr_init2(pi, bits);
  mpfr_mul(pi, c.pi, C, MPFR_RNDN);
  mpfr_clear(c.pi);
}

#ifdef DEBUG
void copy_status(const char *fname) {
  int pid = fork();
  if(pid < 0) {
    printf("Unable to procreate\n");
    return;
  }
  if(pid == 0) {
    char *orig;
    asprintf(&orig, "/proc/%d/status", getppid());
    execl("/usr/bin/cp", "cp", orig, fname, NULL);
    printf("Failed to copy file %s to %s\n", orig, fname);
    exit(-1);
  }
  wait(NULL);
}
#endif /* DEBUG */

int main(int argc, char *argv[]) {
  // Get the precision from the user.
  int precision = 0;
  
  if(argc > 1) {
    precision = atoi(argv[1]);
  }

  if(precision < 1) {
    printf("Enter the precision: ");
    scanf("%d", &precision);
  }

  if(precision < 1) {
    printf("Need precision higher than 0\n");
    return 1;
  }

  // Calculate pi.
  mpfr_t pi;
  calc_pi(pi, precision);

  // Print pi.
  mpfr_printf("Pi to %d digits is:\n%.*Rf\n", precision, precision - 2, pi);

  // Free the MPFR number.
  mpfr_clear(pi);

#ifdef DEBUG
  copy_status("pi_status.txt");
#endif /* DEBUG */

  return 0;
}
