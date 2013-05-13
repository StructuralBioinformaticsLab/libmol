#include <stdio.h>
#include <check.h>

#include "mol.0.0.6/atom.h"

START_TEST (test_create)
{
	mol_atom atom;
	mol_atom_create(&atom, 1);
	fail_if(! (atom.bondis), "mol_atom_create failed");
}
END_TEST

Suite* atom_suite (void) {
	Suite *suite = suite_create("atom");

	TCase *tcase_basic = tcase_create("basic");
	tcase_add_test(tcase_basic, test_create);

	suite_add_tcase(suite, tcase_basic);

	return suite;
}

int main () {
	int number_failed;
	Suite *suite = atom_suite();
	SRunner *runner = srunner_create(suite);
	srunner_run_all(runner, CK_ENV);
	number_failed = srunner_ntests_failed(runner);
	srunner_free(runner);
	return number_failed;
}
