#include <tomfloat.h>

void draw(mp_float *a)
{
   char buf[8192];
   mp_toradix(&(a->mantissa), buf, 10);
   printf("%s * 2^%ld\n", buf, a->exp);
}

int main(void)
{
  mp_float a, b, c, d, e;
  int err;

  mpf_init_multi(100, &a, &b, &c, &d, &e, NULL);

  mpf_const_d(&a, 1); draw(&a);
  mpf_const_d(&b, 2); draw(&b);
  mpf_const_d(&c, 3); draw(&c);
  mpf_const_d(&d, 4); draw(&d);

  mpf_add(&b, &c, &e); printf("2 + 3            == "); draw(&e);
  mpf_sub(&b, &c, &e); printf("2 - 3            ==");  draw(&e);
  mpf_mul(&b, &c, &e); printf("2 * 3            == "); draw(&e);
  mpf_div(&b, &c, &e); printf("2 / 3            == "); draw(&e);
  mpf_add_d(&b, 3, &e); printf("2 + 3            == "); draw(&e);
  mpf_sub_d(&b, 3, &e); printf("2 - 3            ==");  draw(&e);
  mpf_mul_d(&b, 3, &e); printf("2 * 3            == "); draw(&e);
  mpf_div_d	(&b, 3, &e); printf("2 / 3            == "); draw(&e);
  mpf_const_d(&e, 0); mpf_add_d(&e, 1, &e); printf("0 + 1            == "); draw(&e);
  mpf_const_d(&e, 0); mpf_sub_d(&e, 1, &e); printf("0 - 1            == "); draw(&e);
  printf("\n");
  mpf_invsqrt(&d, &e); printf("1/sqrt(4) == 1/2 == "); draw(&e);
  mpf_invsqrt(&c, &e); printf("1/sqrt(3)        == "); draw(&e);
  mpf_inv(&a, &e);     printf("1/1              == "); draw(&e);
  mpf_inv(&b, &e);     printf("1/2              == "); draw(&e);
  mpf_inv(&c, &e);     printf("1/3              == "); draw(&e);
  mpf_inv(&d, &e);     printf("1/4              == "); draw(&e);
  printf("\n");
  mpf_const_pi(&e);    printf("Pi               == "); draw(&e);
  printf("\n");
  mpf_const_e(&e);     printf("e                == "); draw(&e);
  mpf_exp(&c, &e);     printf("e^3              == "); draw(&e);
  mpf_sqrt(&e, &e);    printf("sqrt(e^3)        == "); draw(&e);
  mpf_sqr(&e, &e);     printf("sqrt(e^3)^2      == "); draw(&e);
  printf("\n");
  mpf_cos(&a, &e);     printf("cos(1)           == "); draw(&e);
  mpf_cos(&b, &e);     printf("cos(2)           == "); draw(&e);
  mpf_cos(&c, &e);     printf("cos(3)           == "); draw(&e);
  mpf_cos(&d, &e);     printf("cos(4)           == "); draw(&e);
  mpf_sin(&a, &e);     printf("sin(1)           == "); draw(&e);
  mpf_sin(&b, &e);     printf("sin(2)           == "); draw(&e);
  mpf_sin(&c, &e);     printf("sin(3)           == "); draw(&e);
  mpf_sin(&d, &e);     printf("sin(4)           == "); draw(&e);
  mpf_tan(&a, &e);     printf("tan(1)           == "); draw(&e);
  mpf_tan(&b, &e);     printf("tan(2)           == "); draw(&e);
  mpf_tan(&c, &e);     printf("tan(3)           == "); draw(&e);
  mpf_tan(&d, &e);     printf("tan(4)           == "); draw(&e);
  mpf_inv(&a, &e);     mpf_atan(&e, &e);    printf("atan(1/1)        == "); draw(&e);
  mpf_inv(&b, &e);     mpf_atan(&e, &e);    printf("atan(1/2)        == "); draw(&e);
  mpf_inv(&c, &e);     mpf_atan(&e, &e);    printf("atan(1/3)        == "); draw(&e);
  mpf_inv(&d, &e);     mpf_atan(&e, &e);    printf("atan(1/4)        == "); draw(&e);
  printf("\n");
  #define lntest(x) if ((err = mpf_const_ln_d(&e, x)) != MP_OKAY) { printf("Failed ln(%3d), %d\n", x, err); } else { printf("ln(%3d)          == ", x); draw(&e); };
  lntest(0);
  lntest(1);
  lntest(2);
  lntest(4);
  lntest(8);
  lntest(17);
  lntest(1000);
  lntest(100000);
  lntest(250000);
  return 0; 
}


   
