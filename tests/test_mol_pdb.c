#include <stdio.h>
#include <check.h>
#include <errno.h>

#include "mol.0.0.6/pdb.h"

START_TEST (test_read_small01)
{
	errno = 0;
	mol_atom_group *ag = read_pdb_nopar("small01.pdb");
	ck_assert(ag != NULL);
	ck_assert_int_eq(ag->natoms, 35);
	mol_atom_group_destroy(ag);
}
END_TEST

START_TEST (test_read_small02)
{
	errno = 0;
	mol_atom_group *ag = read_pdb_nopar("small02.pdb");
	ck_assert(ag != NULL);
	ck_assert_int_eq(ag->natoms, 1);
	ck_assert(ag->atoms[0].X == 1142.420);
	ck_assert(ag->atoms[0].Y == 1119.291);
	ck_assert(ag->atoms[0].Z ==  -57.094);
}
END_TEST

Suite *pdb_suite (void) {
	Suite *suite = suite_create("pdb");
	TCase *tcase_basic = tcase_create("basic");
	tcase_add_test(tcase_basic, test_read_small01);
	tcase_add_test(tcase_basic, test_read_small02);

	suite_add_tcase(suite, tcase_basic);

	return suite;
}

int main () {
	Suite *suite = pdb_suite();
	SRunner *runner = srunner_create(suite);
	//srunner_set_fork_status(runner, CK_NOFORK);
	srunner_run_all(runner, CK_ENV);

	int number_failed = srunner_ntests_failed(runner);
	srunner_free(runner);
	return number_failed;
}
