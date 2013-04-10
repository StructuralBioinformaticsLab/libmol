#include <stdio.h>
#include <check.h>

#include "mol.0.0.6/atom.h"

START_TEST (test_create)
{
	mol_atom atom;
	mol_atom_create(&atom, 1);
	fail_if(! (&atom), "atom is should never be null");

}
END_TEST

Suite* atom_suite (void) {
	Suite *suite = suite_create("atom");
	TCase *tcase = tcase_create("case_simple");
	tcase_add_test(tcase, test_create);
	suite_add_tcase(suite, tcase);
	return suite;
}

int main (int argc, char *argv[]) {
	int number_failed;
	Suite *suite = atom_suite();
	SRunner *runner = srunner_create(suite);
	srunner_run_all(runner, CK_NORMAL);
	number_failed = srunner_ntests_failed(runner);
	srunner_free(runner);
	return number_failed;
}
