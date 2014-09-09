#include <stdlib.h>
#include <stdio.h>
#include <check.h>
#include <errno.h>
#include <math.h>

#include "mol.0.0.6/benergy.h"
#include "mol.0.0.6/pdb.h"
#include "mol.0.0.6/icharmm.h"

struct atomgrp *test_ag;
const double delta = 0.000001;
const double tolerance = 0.1;

void check_grads(struct atomgrp *ag, double d, void (*efun)(struct atomgrp *, double*))
{
        int n=ag->natoms, i;
        double en, en1, t;
        double *fs=malloc(3*n*sizeof(double));
//en0
        en=0;
        (*efun)(ag, &en);

        for(i=0; i<n; i++)
        {
//x
                en1=0;
                t=ag->atoms[i].X;
                ag->atoms[i].X=d+t;
                (*efun)(ag, &en1);
                ag->atoms[i].X=t;
                fs[3*i]=(en-en1)/d;
//y
                en1=0;
                t=ag->atoms[i].Y;
                ag->atoms[i].Y=d+t;
                (*efun)(ag, &en1);
                ag->atoms[i].Y=t;
                fs[3*i+1]=(en-en1)/d;
//z
                en1=0;
                t=ag->atoms[i].Z;
                ag->atoms[i].Z=d+t;
                (*efun)(ag, &en1);
                ag->atoms[i].Z=t;
                fs[3*i+2]=(en-en1)/d;
        }
        en=0;
        zero_grads(ag);
        (*efun)(ag, &en);

        char msg[256];
        for(i=0; i<n; i++)
        {
          sprintf(msg,
                  "\n(atom: %d) calc: (%lf, %lf, %lf): numerical: (%lf, %lf, %lf)\n",
                  i,
                  ag->atoms[i].GX, ag->atoms[i].GY, ag->atoms[i].GZ,
                  fs[3*i],fs[3*i+1],fs[3*i+2]);
          ck_assert_msg(fabs(ag->atoms[i].GX - fs[3*i]) < tolerance, msg);
          ck_assert_msg(fabs(ag->atoms[i].GY - fs[3*i+1]) < tolerance, msg);
          ck_assert_msg(fabs(ag->atoms[i].GZ - fs[3*i+2]) < tolerance, msg);

        }
        free(fs);
}

void setup(void)
{
	test_ag = read_pdb_nopar("1rei_nmin.pdb");
	read_ff_charmm("1rei_nmin.psf", "parm.prm", "pdbamino.rtf", test_ag);
	fixed_init(test_ag);
	fixed_update(test_ag, 0, NULL); // Make nothing fixed
}

void teardown(void)
{

}

// Test cases
START_TEST(test_beng)
{
	zero_grads(test_ag);
	//check_grads(test_ag, delta, beng);
}
END_TEST

START_TEST(test_aeng)
{
	zero_grads(test_ag);
	//check_grads(test_ag, delta, aeng);
}
END_TEST

START_TEST(test_ieng)
{
	zero_grads(test_ag);
	check_grads(test_ag, delta, ieng);
}
END_TEST

START_TEST(test_teng)
{
	check_grads(test_ag, delta, teng);
}
END_TEST

Suite *benergy_suite(void)
{
	Suite *suite = suite_create("benergy");

	// Add test cases here
	// Each test case can call multiple test functions
	// Too add a test case, call tcase_add_test
	// The first argument is the TCase struct, the second is the
	//  test function name.
	TCase *tcase = tcase_create("test");
	tcase_set_timeout(tcase, 20);
        tcase_add_checked_fixture(tcase, setup, teardown);
	tcase_add_test(tcase, test_beng);
	tcase_add_test(tcase, test_aeng);
	tcase_add_test(tcase, test_ieng);
	tcase_add_test(tcase, test_teng);

	suite_add_tcase(suite, tcase);

	return suite;
}

int main(void)
{
	Suite *suite = benergy_suite();
	SRunner *runner = srunner_create(suite);
	srunner_run_all(runner, CK_ENV);

	int number_failed = srunner_ntests_failed(runner);
	srunner_free(runner);
	return number_failed;
}
