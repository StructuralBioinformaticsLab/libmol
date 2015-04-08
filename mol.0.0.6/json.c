/*
Copyright (c) 2013, Acpharis
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
- Neither the name of the author nor the names of its contributors may be used
  to endorse or promote products derived from this software without specific
  prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#ifndef _NO_JANSSON_
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <jansson.h>
#include _MOL_INCLUDE_

struct atomgrp* read_json_ag(const char *json_file)
{
        struct atomgrp* ag = (struct atomgrp*) _mol_calloc(1, sizeof (struct atomgrp));
	json_error_t json_file_error;
	json_t *base = json_load_file(json_file, 0, &json_file_error); 

	if (!base) {
		fprintf(stderr, "error reading json file %s on line %d column %d: %s\n", json_file, json_file_error.line, json_file_error.column, json_file_error.text);
	}

	if (!json_is_object(base)) {
		fprintf(stderr, "json file not an object %s\n", json_file);
	}

	json_t *atoms, *bonds, *angles, *torsions, *impropers;
	atoms = json_object_get(base, "atoms");
	if (!json_is_array(atoms)) {
		fprintf(stderr, "json atoms are not an array %s\n", json_file);
	}
	size_t natoms = json_array_size(atoms);
	ag->natoms = natoms;
        ag->atoms = (struct atom*) _mol_calloc(ag->natoms, sizeof(struct atom));

	ag->num_atom_types = 0;

	char *prev_segment = _mol_calloc(1, sizeof(char));
	char *prev_residue = _mol_calloc(1, sizeof(char));
	int prev_residue_seq = -107;
	for (size_t i=0; i < natoms; i++) {
		json_t *atom = json_array_get(atoms, i);
		if (!json_is_object(atom)) {
			fprintf(stderr, "Atom %zd not an object in json file %s\n", i, json_file);
		}
		json_t *ace_volume, *ftype_index, *ftype_name, *eps03;
		json_t *name, *radius03, *eps, *acp_type, *residue_name;
		json_t *charge, *radius, *element;
		json_t *x, *y, *z;
		json_t *yeti_type, *sybyl_type;
		json_t *backbone, *hb_acceptor, *hb_donor, *hb_weight;
		json_t *segment, *residue;

		segment = json_object_get(atom, "segment");
		residue = json_object_get(atom, "residue");
		if ((segment != NULL) && (residue != NULL)) {
			if (!json_is_string(segment)) {
				fprintf(stderr, "json segment is not string for atom %zd in json_file %s\n", i, json_file);
			}
			if (!json_is_string(residue)) {
				fprintf(stderr, "json residue is not string for atom %zd in json_file %s\n", i, json_file);
			}

			const char *cur_segment = json_string_value(segment);
			const char *cur_residue = json_string_value(residue);

			if (strcmp(cur_segment, prev_segment) != 0) {
				prev_residue_seq += 100;
				free(prev_segment);
				prev_segment = strdup(cur_segment);
			}
			if (strcmp(cur_residue, prev_residue) != 0) {
				int cur_residue_int = atoi(cur_residue);
				int prev_residue_int = atoi(prev_residue);
				if ( (cur_residue_int-prev_residue_int) > 1) {
					prev_residue_seq += (cur_residue_int-prev_residue_int);
				} else {
					prev_residue_seq += 1;
				}
				free(prev_residue);
				prev_residue = strdup(cur_residue);
			}

			ag->atoms[i].comb_res_seq = prev_residue_seq;
		} else {
			ag->atoms[i].comb_res_seq = prev_residue_seq;
		}

		ace_volume = json_object_get(atom, "ace_volume");
		if (!json_is_real(ace_volume)) {
			fprintf(stderr, "json ace volume is not floating point for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].acevolume = json_real_value(ace_volume);

		ftype_index = json_object_get(atom, "ftype_index");
		if (!json_is_integer(ftype_index)) {
			fprintf(stderr, "json ftype index is not integer for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].atom_ftypen = json_integer_value(ftype_index);
		if (ag->atoms[i].atom_ftypen > ag->num_atom_types) {
			ag->num_atom_types = ag->atoms[i].atom_ftypen;
		}

		ftype_name = json_object_get(atom, "ftype_name");
		if (!json_is_string(ftype_name)) {
			fprintf(stderr, "json ftype name is not string for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].ftype_name = strdup(json_string_value(ftype_name));

		element = json_object_get(atom, "element");
		if (!json_is_string(element)) {
			fprintf(stderr, "json element name is not string for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].element = strdup(json_string_value(element));

		eps = json_object_get(atom, "eps");
		if (!json_is_real(eps)) {
			fprintf(stderr, "json eps is not floating point for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].eps = sqrt(-json_real_value(eps));

		eps03 = json_object_get(atom, "eps03");
		if (!json_is_real(eps03)) {
			fprintf(stderr, "json eps03 is not floating point for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].eps03 = sqrt(-json_real_value(eps03));

		radius = json_object_get(atom, "radius");
		if (!json_is_real(radius)) {
			fprintf(stderr, "json radius is not floating point for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].rminh = json_real_value(radius);

		radius03 = json_object_get(atom, "radius03");
		if (!json_is_real(radius03)) {
			fprintf(stderr, "json radius03 is not floating point for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].rminh03 = json_real_value(radius03);

		charge = json_object_get(atom, "charge");
		if (!json_is_real(radius03)) {
			fprintf(stderr, "json charge is not floating point for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].chrg = json_real_value(charge);

		name = json_object_get(atom, "name");
		if (!json_is_string(name)) {
			fprintf(stderr, "json name is not string for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].name = strdup(json_string_value(name));

		residue_name = json_object_get(atom, "residue_name");
		if (residue_name != NULL) {
			if (!json_is_string(residue_name)) {
				fprintf(stderr, "json residue_name is not string for atom %zd in json_file %s\n", i, json_file);
			}
		}
		ag->atoms[i].residue_name = strdup(json_string_value(residue_name));

		x = json_object_get(atom, "x");
		if (!json_is_real(x)) {
			fprintf(stderr, "json coordinate x is not floating point for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].X = json_real_value(x);

		y = json_object_get(atom, "y");
		if (!json_is_real(y)) {
			fprintf(stderr, "json coordinate y is not floating point for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].Y = json_real_value(y);

		z = json_object_get(atom, "z");
		if (!json_is_real(z)) {
			fprintf(stderr, "json coordinate z is not floating point for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].Z = json_real_value(z);

		acp_type = json_object_get(atom, "acp_type");
		if (acp_type != NULL) {
			if (!json_is_integer(acp_type)) {
				fprintf(stderr, "json acp_type index is not integer for atom %zd in json_file %s\n", i, json_file);
			}
			ag->atoms[i].acp_type = json_integer_value(acp_type);
		} else {
			ag->atoms[i].acp_type = -1;
		}

		yeti_type = json_object_get(atom, "yeti_type");
		if (yeti_type != NULL) {
			const char *yeti_type_string;
			if (!json_is_string(yeti_type)) {
				fprintf(stderr, "json yeti_type is not string for atom %zd in json_file %s\n", i, json_file);
			}
			yeti_type_string = json_string_value(yeti_type);
			if (strcmp("carbonyl", yeti_type_string) == 0) {
				ag->atoms[i].yeti_type = MOL_YETI_CARBONYL;
			} else if (strcmp("hydroxyl", yeti_type_string) == 0) {
				ag->atoms[i].yeti_type = MOL_YETI_HYDROXYL;
			} else if (strcmp("sulfonamide", yeti_type_string) == 0) {
				ag->atoms[i].yeti_type = MOL_YETI_SULFONAMIDE;
			} else if (strcmp("N5_aromatic", yeti_type_string) == 0) {
				ag->atoms[i].yeti_type = MOL_YETI_N5_AROMATIC;
			} else if (strcmp("N6_aromatic", yeti_type_string) == 0) {
				ag->atoms[i].yeti_type = MOL_YETI_N6_AROMATIC;
			} else {
				fprintf(stderr, "unknown json yeti_type %s for atom %zd in json_file %s\n", yeti_type_string, i, json_file);
				ag->atoms[i].yeti_type = MOL_YETI_NONE;
			}
		} else {
			ag->atoms[i].yeti_type = MOL_YETI_NONE;
		}

		sybyl_type = json_object_get(atom, "sybyl_type");
		if (sybyl_type != NULL) {
			const char *sybyl_type_string;
			if (!json_is_string(sybyl_type)) {
				fprintf(stderr, "json sybyl_type is not string for atom %zd in json_file %s\n", i, json_file);
			}
			sybyl_type_string = json_string_value(sybyl_type);
			ag->atoms[i].hybridization = mol_hydridization_from_sybyl(sybyl_type_string);
		} else {
			ag->atoms[i].hybridization = UNKNOWN_HYBRID;
		}

		backbone = json_object_get(atom, "backbone");
		if (backbone != NULL) {
			if (!json_is_boolean(backbone)) {
				fprintf(stderr, "json backbone is not boolean for atom %zd in json_file %s\n", i, json_file);
			}
			ag->atoms[i].backbone = json_is_true(backbone);
		} else {
			ag->atoms[i].backbone = false;
		}

		ag->atoms[i].hprop = 0;
		hb_acceptor = json_object_get(atom, "hb_acceptor");
		if (hb_acceptor != NULL) {
			if (!json_is_boolean(hb_acceptor)) {
				fprintf(stderr, "json hb_acceptor is not boolean for atom %zd in json_file %s\n", i, json_file);
			}
			if (json_is_true(hb_acceptor)) {
				ag->atoms[i].hprop |= HBOND_ACCEPTOR;
			}
		}
		hb_donor = json_object_get(atom, "hb_donor");
		if (hb_donor != NULL) {
			if (!json_is_boolean(hb_donor)) {
				fprintf(stderr, "json hb_donor is not boolean for atom %zd in json_file %s\n", i, json_file);
			}
			if (json_is_true(hb_donor)) {
				ag->atoms[i].hprop |= DONATABLE_HYDROGEN;
			}
		}
		ag->atoms[i].hbond_weight = 1.0;
		hb_weight = json_object_get(atom, "hb_weight");
		if (hb_weight != NULL) {
			if (!json_is_real(hb_weight)) {
				fprintf(stderr, "json hb_weight is not floating point for atom %zd in json_file %s\n", i, json_file);
			}
			ag->atoms[i].hbond_weight = json_real_value(hb_weight);
		}

		ag->atoms[i].nbonds = 0;
		ag->atoms[i].nangs = 0;
		ag->atoms[i].ntors = 0;
		ag->atoms[i].nimps = 0;
		ag->atoms[i].fixed = 0;
		ag->atoms[i].sa = -1;
		ag->atoms[i].base = -1;
		ag->atoms[i].base2 = -1;
	}
	free(prev_segment);
	free(prev_residue);

	bonds = json_object_get(base, "bonds");
	if (!json_is_array(bonds)) {
		fprintf(stderr, "json bonds are not an array %s\n", json_file);
	}
	size_t nbonds = json_array_size(bonds);
	ag->nbonds = nbonds;
	ag->bonds = _mol_calloc(nbonds, sizeof(struct atombond));
	for (size_t i=0; i < nbonds; i++) {
		json_t *bond = json_array_get(bonds, i);
		if (!json_is_object(bond)) {
			fprintf(stderr, "Bond %zd not an object in json file %s\n", i, json_file);
		}
		json_t *length, *atom1, *atom2, *spring_constant, *sdf_type;

		atom1 = json_object_get(bond, "atom1");
		if (!json_is_integer(atom1)) {
			fprintf(stderr, "json atom1 is not integer for bond %zd in json_file %s\n", i, json_file);
		}
		int i1 = json_integer_value(atom1) - 1;
		ag->bonds[i].a0 = &(ag->atoms[i1]);
		(ag->atoms[i1].nbonds)++;

		atom2 = json_object_get(bond, "atom2");
		if (!json_is_integer(atom2)) {
			fprintf(stderr, "json atom2 is not integer for bond %zd in json_file %s\n", i, json_file);
		}
		int i2 = json_integer_value(atom2) - 1;
		ag->bonds[i].a1 = &(ag->atoms[i2]);
		(ag->atoms[i2].nbonds)++;

		length = json_object_get(bond, "length");
		if (!json_is_real(length)) {
			fprintf(stderr, "json length is not floating point for bond %zd in json_file %s\n", i, json_file);
		}
		ag->bonds[i].l0 = json_real_value(length);

		spring_constant = json_object_get(bond, "spring_constant");
		if (!json_is_real(spring_constant)) {
			fprintf(stderr, "json spring_constant is not floating point for bond %zd in json_file %s\n", i, json_file);
		}
		ag->bonds[i].k = json_real_value(spring_constant);

		sdf_type = json_object_get(bond, "sdf_type");
		if (sdf_type != NULL) {
			if (!json_is_integer(sdf_type)) {
				fprintf(stderr, "json sdf_type is not integer for bond %zd in json_file %s\n", i, json_file);
			}
			ag->bonds[i].sdf_type = json_integer_value(sdf_type);
		} else {
			ag->bonds[i].sdf_type = 0;
		}
	}

	angles = json_object_get(base, "angles");
	if (!json_is_array(angles)) {
		fprintf(stderr, "json angles are not an array %s\n", json_file);
	}
	size_t nangles = json_array_size(angles);
	ag->nangs = nangles;
	ag->angs = _mol_calloc(nangles, sizeof(struct atomangle));
	for (size_t i=0; i < nangles; i++) {
		json_t *angle = json_array_get(angles, i);
		if (!json_is_object(angle)) {
			fprintf(stderr, "Angle %zd not an object in json file %s\n", i, json_file);
		}
		json_t *theta, *atom1, *atom2, *atom3, *spring_constant;

		atom1 = json_object_get(angle, "atom1");
		if (!json_is_integer(atom1)) {
			fprintf(stderr, "json atom1 is not integer for angle %zd in json_file %s\n", i, json_file);
		}
		int i1 = json_integer_value(atom1) - 1;
		ag->angs[i].a0 = &(ag->atoms[i1]);
		(ag->atoms[i1].nangs)++;

		atom2 = json_object_get(angle, "atom2");
		if (!json_is_integer(atom2)) {
			fprintf(stderr, "json atom2 is not integer for angle %zd in json_file %s\n", i, json_file);
		}
		int i2 = json_integer_value(atom2) - 1;
		ag->angs[i].a1 = &(ag->atoms[i2]);
		(ag->atoms[i2].nangs)++;

		atom3 = json_object_get(angle, "atom3");
		if (!json_is_integer(atom3)) {
			fprintf(stderr, "json atom3 is not integer for angle %zd in json_file %s\n", i, json_file);
		}
		int i3 = json_integer_value(atom3) - 1;
		ag->angs[i].a2 = &(ag->atoms[i3]);
		(ag->atoms[i3].nangs)++;

		theta = json_object_get(angle, "theta");
		if (!json_is_real(theta)) {
			fprintf(stderr, "json theta is not floating point for angle %zd in json_file %s\n", i, json_file);
		}
		ag->angs[i].th0 = json_real_value(theta);

		spring_constant = json_object_get(angle, "spring_constant");
		if (!json_is_real(spring_constant)) {
			fprintf(stderr, "json spring_constant is not floating point for angle %zd in json_file %s\n", i, json_file);
		}
		ag->angs[i].k = json_real_value(spring_constant);
	}

	torsions = json_object_get(base, "torsions");
	if (!json_is_array(torsions)) {
		fprintf(stderr, "json torsions are not an array %s\n", json_file);
	}
	size_t ntorsions = json_array_size(torsions);
	ag->ntors = ntorsions;
	ag->tors = _mol_calloc(ntorsions, sizeof(struct atomtorsion));
	for (size_t i=0; i < ntorsions; i++) {
		json_t *torsion = json_array_get(torsions, i);
		if (!json_is_object(torsion)) {
			fprintf(stderr, "Torsion %zd not an object in json file %s\n", i, json_file);
		}
		json_t *atom1, *atom2, *atom3, *atom4, *minima, *delta_constant, *spring_constant;

		atom1 = json_object_get(torsion, "atom1");
		if (!json_is_integer(atom1)) {
			fprintf(stderr, "json atom1 is not integer for torsion %zd in json_file %s\n", i, json_file);
		}
		int i1 = json_integer_value(atom1) - 1;
		ag->tors[i].a0 = &(ag->atoms[i1]);
		(ag->atoms[i1].ntors)++;

		atom2 = json_object_get(torsion, "atom2");
		if (!json_is_integer(atom2)) {
			fprintf(stderr, "json atom2 is not integer for torsion %zd in json_file %s\n", i, json_file);
		}
		int i2 = json_integer_value(atom2) - 1;
		ag->tors[i].a1 = &(ag->atoms[i2]);
		(ag->atoms[i2].ntors)++;

		atom3 = json_object_get(torsion, "atom3");
		if (!json_is_integer(atom3)) {
			fprintf(stderr, "json atom3 is not integer for torsion %zd in json_file %s\n", i, json_file);
		}
		int i3 = json_integer_value(atom3) - 1;
		ag->tors[i].a2 = &(ag->atoms[i3]);
		(ag->atoms[i3].ntors)++;

		atom4 = json_object_get(torsion, "atom4");
		if (!json_is_integer(atom4)) {
			fprintf(stderr, "json atom4 is not integer for torsion %zd in json_file %s\n", i, json_file);
		}
		int i4 = json_integer_value(atom4) - 1;
		ag->tors[i].a3 = &(ag->atoms[i4]);
		(ag->atoms[i4].ntors)++;

		minima = json_object_get(torsion, "minima");
		if (!json_is_integer(minima)) {
			fprintf(stderr, "json minima is not integer for torsion %zd in json_file %s\n", i, json_file);
		}
		ag->tors[i].n = json_integer_value(minima);

		delta_constant = json_object_get(torsion, "delta_constant");
		if (!json_is_real(delta_constant)) {
			fprintf(stderr, "json delta_constant is not floating point for torsion %zd in json_file %s\n", i, json_file);
		}
		ag->tors[i].d = json_real_value(delta_constant);

		spring_constant = json_object_get(torsion, "spring_constant");
		if (!json_is_real(spring_constant)) {
			fprintf(stderr, "json spring_constant is not floating point for torsion %zd in json_file %s\n", i, json_file);
		}
		ag->tors[i].k = json_real_value(spring_constant);
	}

	impropers = json_object_get(base, "impropers");
	if (!json_is_array(impropers)) {
		fprintf(stderr, "json impropers are not an array %s\n", json_file);
	}
	size_t nimpropers = json_array_size(impropers);
	ag->nimps = nimpropers;
	ag->imps = _mol_calloc(nimpropers, sizeof(struct atomimproper));
	for (size_t i=0; i < nimpropers; i++) {
		json_t *improper = json_array_get(impropers, i);
		if (!json_is_object(improper)) {
			fprintf(stderr, "Improper %zd not an object in json file %s\n", i, json_file);
		}
		json_t *atom1, *atom2, *atom3, *atom4, *phi, *spring_constant;

		atom1 = json_object_get(improper, "atom1");
		if (!json_is_integer(atom1)) {
			fprintf(stderr, "json atom1 is not integer for improper %zd in json_file %s\n", i, json_file);
		}
		int i1 = json_integer_value(atom1) - 1;
		ag->imps[i].a0 = &(ag->atoms[i1]);
		(ag->atoms[i1].nimps)++;

		atom2 = json_object_get(improper, "atom2");
		if (!json_is_integer(atom2)) {
			fprintf(stderr, "json atom2 is not integer for improper %zd in json_file %s\n", i, json_file);
		}
		int i2 = json_integer_value(atom2) - 1;
		ag->imps[i].a1 = &(ag->atoms[i2]);
		(ag->atoms[i2].nimps)++;

		atom3 = json_object_get(improper, "atom3");
		if (!json_is_integer(atom3)) {
			fprintf(stderr, "json atom3 is not integer for improper %zd in json_file %s\n", i, json_file);
		}
		int i3 = json_integer_value(atom3) - 1;
		ag->imps[i].a2 = &(ag->atoms[i3]);
		(ag->atoms[i3].nimps)++;

		atom4 = json_object_get(improper, "atom4");
		if (!json_is_integer(atom4)) {
			fprintf(stderr, "json atom4 is not integer for improper %zd in json_file %s\n", i, json_file);
		}
		int i4 = json_integer_value(atom4) - 1;
		ag->imps[i].a3 = &(ag->atoms[i4]);
		(ag->atoms[i4].nimps)++;

		phi = json_object_get(improper, "phi");
		if (!json_is_real(phi)) {
			fprintf(stderr, "json phi is not floating point for improper %zd in json_file %s\n", i, json_file);
		}
		ag->imps[i].psi0 = json_real_value(phi);

		spring_constant = json_object_get(improper, "spring_constant");
		if (!json_is_real(spring_constant)) {
			fprintf(stderr, "json spring_constant is not floating point for improper %zd in json_file %s\n", i, json_file);
		}
		ag->imps[i].k = json_real_value(spring_constant);
	}



	json_decref(base);

//allocate atom arrays of pointers to parameters
	for (size_t i = 0; i < natoms; i++) {
		int i1 = ag->atoms[i].nbonds;
		ag->atoms[i].bonds =
		    _mol_calloc(i1, sizeof(struct atombond *));
		ag->atoms[i].nbonds = 0;
		i1 = ag->atoms[i].nangs;
		ag->atoms[i].angs =
		    _mol_calloc(i1, sizeof(struct atomangle *));
		ag->atoms[i].nangs = 0;
		i1 = ag->atoms[i].ntors;
		ag->atoms[i].tors =
		    _mol_calloc(i1, sizeof(struct atomtorsion *));
		ag->atoms[i].ntors = 0;
		i1 = ag->atoms[i].nimps;
		ag->atoms[i].imps =
		    _mol_calloc(i1, sizeof(struct atomimproper *));
		ag->atoms[i].nimps = 0;
	}
	struct atom *atm;
//fill bonds
	for (int i = 0; i < ag->nbonds; i++) {
		atm = ag->bonds[i].a0;
		atm->bonds[(atm->nbonds)++] = &(ag->bonds[i]);
		atm = ag->bonds[i].a1;
		atm->bonds[(atm->nbonds)++] = &(ag->bonds[i]);
	}
//fill angles
	for (int i = 0; i < ag->nangs; i++) {
		atm = ag->angs[i].a0;
		atm->angs[(atm->nangs)++] = &(ag->angs[i]);
		atm = ag->angs[i].a1;
		atm->angs[(atm->nangs)++] = &(ag->angs[i]);
		atm = ag->angs[i].a2;
		atm->angs[(atm->nangs)++] = &(ag->angs[i]);
	}
//fill torsions
	for (int i = 0; i < ag->ntors; i++) {
		atm = ag->tors[i].a0;
		atm->tors[(atm->ntors)++] = &(ag->tors[i]);
		atm = ag->tors[i].a1;
		atm->tors[(atm->ntors)++] = &(ag->tors[i]);
		atm = ag->tors[i].a2;
		atm->tors[(atm->ntors)++] = &(ag->tors[i]);
		atm = ag->tors[i].a3;
		atm->tors[(atm->ntors)++] = &(ag->tors[i]);
	}
//fill impropers
	for (int i = 0; i < ag->nimps; i++) {
		atm = ag->imps[i].a0;
		atm->imps[(atm->nimps)++] = &(ag->imps[i]);
		atm = ag->imps[i].a1;
		atm->imps[(atm->nimps)++] = &(ag->imps[i]);
		atm = ag->imps[i].a2;
		atm->imps[(atm->nimps)++] = &(ag->imps[i]);
		atm = ag->imps[i].a3;
		atm->imps[(atm->nimps)++] = &(ag->imps[i]);
	}
//atom indices in the group
	fill_ingrp(ag);

	ag->is_psf_read = true;

	// copy vals from deprecated to new data structures
	int atomsi;
	for (atomsi = 0; atomsi < ag->natoms; atomsi++) {
		_mol_atom_create_bond_indices(&ag->atoms[atomsi],
					      ag->atoms[atomsi].nbonds);
	}
	_mol_atom_group_copy_from_deprecated(ag);
	return ag;
}


void read_ff_json(const char *json_file, struct atomgrp *ag)
{
	json_error_t json_file_error;
	json_t *base = json_load_file(json_file, 0, &json_file_error); 

	if (!base) {
		fprintf(stderr, "error reading json file %s on line %d column %d: %s\n", json_file, json_file_error.line, json_file_error.column, json_file_error.text);
	}

	if (!json_is_object(base)) {
		fprintf(stderr, "json file not an object %s\n", json_file);
	}

	json_t *atoms, *bonds, *angles, *torsions, *impropers;
	atoms = json_object_get(base, "atoms");
	if (!json_is_array(atoms)) {
		fprintf(stderr, "json atoms are not an array %s\n", json_file);
	}
	size_t natoms = json_array_size(atoms);
	if (natoms != (size_t) ag->natoms) {
		fprintf(stderr, "json file has a different number of atoms %zd vs. %d : %s\n", natoms, ag->natoms, json_file);
	}

	ag->num_atom_types = 0;
	for (size_t i=0; i < natoms; i++) {
		json_t *atom = json_array_get(atoms, i);
		if (!json_is_object(atom)) {
			fprintf(stderr, "Atom %zd not an object in json file %s\n", i, json_file);
		}
		json_t *ace_volume, *ftype_index, *ftype_name, *eps03;
		json_t *name, *radius03, *eps;
		json_t *charge, *radius, *element;

		ace_volume = json_object_get(atom, "ace_volume");
		if (!json_is_real(ace_volume)) {
			fprintf(stderr, "json ace volume is not floating point for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].acevolume = json_real_value(ace_volume);

		ftype_index = json_object_get(atom, "ftype_index");
		if (!json_is_integer(ftype_index)) {
			fprintf(stderr, "json ftype index is not integer for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].atom_ftypen = json_integer_value(ftype_index);
		if (ag->atoms[i].atom_ftypen > ag->num_atom_types) {
			ag->num_atom_types = ag->atoms[i].atom_ftypen;
		}

		ftype_name = json_object_get(atom, "ftype_name");
		if (!json_is_string(ftype_name)) {
			fprintf(stderr, "json ftype name is not string for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].ftype_name = strdup(json_string_value(ftype_name));

		element = json_object_get(atom, "element");
		if (!json_is_string(element)) {
			fprintf(stderr, "json element name is not string for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].element = strdup(json_string_value(element));

		eps = json_object_get(atom, "eps");
		if (!json_is_real(eps)) {
			fprintf(stderr, "json eps is not floating point for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].eps = sqrt(-json_real_value(eps));

		eps03 = json_object_get(atom, "eps03");
		if (!json_is_real(eps03)) {
			fprintf(stderr, "json eps03 is not floating point for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].eps03 = sqrt(-json_real_value(eps03));

		radius = json_object_get(atom, "radius");
		if (!json_is_real(radius)) {
			fprintf(stderr, "json radius is not floating point for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].rminh = json_real_value(radius);

		radius03 = json_object_get(atom, "radius03");
		if (!json_is_real(radius03)) {
			fprintf(stderr, "json radius03 is not floating point for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].rminh03 = json_real_value(radius03);

		charge = json_object_get(atom, "charge");
		if (!json_is_real(radius03)) {
			fprintf(stderr, "json charge is not floating point for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].chrg = json_real_value(charge);

		name = json_object_get(atom, "name");
		if (!json_is_string(name)) {
			fprintf(stderr, "json name is not string for atom %zd in json_file %s\n", i, json_file);
		}
		ag->atoms[i].name = strdup(json_string_value(name));

		ag->atoms[i].nbonds = 0;
		ag->atoms[i].nangs = 0;
		ag->atoms[i].ntors = 0;
		ag->atoms[i].nimps = 0;
		ag->atoms[i].fixed = 0;
	}

	bonds = json_object_get(base, "bonds");
	if (!json_is_array(bonds)) {
		fprintf(stderr, "json bonds are not an array %s\n", json_file);
	}
	size_t nbonds = json_array_size(bonds);
	ag->nbonds = nbonds;
	ag->bonds = _mol_calloc(nbonds, sizeof(struct atombond));
	for (size_t i=0; i < nbonds; i++) {
		json_t *bond = json_array_get(bonds, i);
		if (!json_is_object(bond)) {
			fprintf(stderr, "Bond %zd not an object in json file %s\n", i, json_file);
		}
		json_t *length, *atom1, *atom2, *spring_constant, *sdf_type;

		atom1 = json_object_get(bond, "atom1");
		if (!json_is_integer(atom1)) {
			fprintf(stderr, "json atom1 is not integer for bond %zd in json_file %s\n", i, json_file);
		}
		int i1 = json_integer_value(atom1) - 1;
		ag->bonds[i].a0 = &(ag->atoms[i1]);
		(ag->atoms[i1].nbonds)++;

		atom2 = json_object_get(bond, "atom2");
		if (!json_is_integer(atom2)) {
			fprintf(stderr, "json atom2 is not integer for bond %zd in json_file %s\n", i, json_file);
		}
		int i2 = json_integer_value(atom2) - 1;
		ag->bonds[i].a1 = &(ag->atoms[i2]);
		(ag->atoms[i2].nbonds)++;

		length = json_object_get(bond, "length");
		if (!json_is_real(length)) {
			fprintf(stderr, "json length is not floating point for bond %zd in json_file %s\n", i, json_file);
		}
		ag->bonds[i].l0 = json_real_value(length);

		spring_constant = json_object_get(bond, "spring_constant");
		if (!json_is_real(spring_constant)) {
			fprintf(stderr, "json spring_constant is not floating point for bond %zd in json_file %s\n", i, json_file);
		}
		ag->bonds[i].k = json_real_value(spring_constant);

		sdf_type = json_object_get(bond, "sdf_type");
		if (!json_is_integer(sdf_type)) {
			fprintf(stderr, "json sdf_type is not integer for bond %zd in json_file %s\n", i, json_file);
		}
		ag->bonds[i].sdf_type = json_integer_value(sdf_type);
	}

	angles = json_object_get(base, "angles");
	if (!json_is_array(angles)) {
		fprintf(stderr, "json angles are not an array %s\n", json_file);
	}
	size_t nangles = json_array_size(angles);
	ag->nangs = nangles;
	ag->angs = _mol_calloc(nangles, sizeof(struct atomangle));
	for (size_t i=0; i < nangles; i++) {
		json_t *angle = json_array_get(angles, i);
		if (!json_is_object(angle)) {
			fprintf(stderr, "Angle %zd not an object in json file %s\n", i, json_file);
		}
		json_t *theta, *atom1, *atom2, *atom3, *spring_constant;

		atom1 = json_object_get(angle, "atom1");
		if (!json_is_integer(atom1)) {
			fprintf(stderr, "json atom1 is not integer for angle %zd in json_file %s\n", i, json_file);
		}
		int i1 = json_integer_value(atom1) - 1;
		ag->angs[i].a0 = &(ag->atoms[i1]);
		(ag->atoms[i1].nangs)++;

		atom2 = json_object_get(angle, "atom2");
		if (!json_is_integer(atom2)) {
			fprintf(stderr, "json atom2 is not integer for angle %zd in json_file %s\n", i, json_file);
		}
		int i2 = json_integer_value(atom2) - 1;
		ag->angs[i].a1 = &(ag->atoms[i2]);
		(ag->atoms[i2].nangs)++;

		atom3 = json_object_get(angle, "atom3");
		if (!json_is_integer(atom3)) {
			fprintf(stderr, "json atom3 is not integer for angle %zd in json_file %s\n", i, json_file);
		}
		int i3 = json_integer_value(atom3) - 1;
		ag->angs[i].a2 = &(ag->atoms[i3]);
		(ag->atoms[i3].nangs)++;

		theta = json_object_get(angle, "theta");
		if (!json_is_real(theta)) {
			fprintf(stderr, "json theta is not floating point for angle %zd in json_file %s\n", i, json_file);
		}
		ag->angs[i].th0 = json_real_value(theta);

		spring_constant = json_object_get(angle, "spring_constant");
		if (!json_is_real(spring_constant)) {
			fprintf(stderr, "json spring_constant is not floating point for angle %zd in json_file %s\n", i, json_file);
		}
		ag->angs[i].k = json_real_value(spring_constant);
	}

	torsions = json_object_get(base, "torsions");
	if (!json_is_array(torsions)) {
		fprintf(stderr, "json torsions are not an array %s\n", json_file);
	}
	size_t ntorsions = json_array_size(torsions);
	ag->ntors = ntorsions;
	ag->tors = _mol_calloc(ntorsions, sizeof(struct atomtorsion));
	for (size_t i=0; i < ntorsions; i++) {
		json_t *torsion = json_array_get(torsions, i);
		if (!json_is_object(torsion)) {
			fprintf(stderr, "Torsion %zd not an object in json file %s\n", i, json_file);
		}
		json_t *atom1, *atom2, *atom3, *atom4, *minima, *delta_constant, *spring_constant;

		atom1 = json_object_get(torsion, "atom1");
		if (!json_is_integer(atom1)) {
			fprintf(stderr, "json atom1 is not integer for torsion %zd in json_file %s\n", i, json_file);
		}
		int i1 = json_integer_value(atom1) - 1;
		ag->tors[i].a0 = &(ag->atoms[i1]);
		(ag->atoms[i1].ntors)++;

		atom2 = json_object_get(torsion, "atom2");
		if (!json_is_integer(atom2)) {
			fprintf(stderr, "json atom2 is not integer for torsion %zd in json_file %s\n", i, json_file);
		}
		int i2 = json_integer_value(atom2) - 1;
		ag->tors[i].a1 = &(ag->atoms[i2]);
		(ag->atoms[i2].ntors)++;

		atom3 = json_object_get(torsion, "atom3");
		if (!json_is_integer(atom3)) {
			fprintf(stderr, "json atom3 is not integer for torsion %zd in json_file %s\n", i, json_file);
		}
		int i3 = json_integer_value(atom3) - 1;
		ag->tors[i].a2 = &(ag->atoms[i3]);
		(ag->atoms[i3].ntors)++;

		atom4 = json_object_get(torsion, "atom4");
		if (!json_is_integer(atom4)) {
			fprintf(stderr, "json atom4 is not integer for torsion %zd in json_file %s\n", i, json_file);
		}
		int i4 = json_integer_value(atom4) - 1;
		ag->tors[i].a3 = &(ag->atoms[i4]);
		(ag->atoms[i4].ntors)++;

		minima = json_object_get(torsion, "minima");
		if (!json_is_integer(minima)) {
			fprintf(stderr, "json minima is not integer for torsion %zd in json_file %s\n", i, json_file);
		}
		ag->tors[i].n = json_integer_value(minima);

		delta_constant = json_object_get(torsion, "delta_constant");
		if (!json_is_real(delta_constant)) {
			fprintf(stderr, "json delta_constant is not floating point for torsion %zd in json_file %s\n", i, json_file);
		}
		ag->tors[i].d = json_real_value(delta_constant);

		spring_constant = json_object_get(torsion, "spring_constant");
		if (!json_is_real(spring_constant)) {
			fprintf(stderr, "json spring_constant is not floating point for torsion %zd in json_file %s\n", i, json_file);
		}
		ag->tors[i].k = json_real_value(spring_constant);
	}

	impropers = json_object_get(base, "impropers");
	if (!json_is_array(impropers)) {
		fprintf(stderr, "json impropers are not an array %s\n", json_file);
	}
	size_t nimpropers = json_array_size(impropers);
	ag->nimps = nimpropers;
	ag->imps = _mol_calloc(nimpropers, sizeof(struct atomimproper));
	for (size_t i=0; i < nimpropers; i++) {
		json_t *improper = json_array_get(impropers, i);
		if (!json_is_object(improper)) {
			fprintf(stderr, "Improper %zd not an object in json file %s\n", i, json_file);
		}
		json_t *atom1, *atom2, *atom3, *atom4, *phi, *spring_constant;

		atom1 = json_object_get(improper, "atom1");
		if (!json_is_integer(atom1)) {
			fprintf(stderr, "json atom1 is not integer for improper %zd in json_file %s\n", i, json_file);
		}
		int i1 = json_integer_value(atom1) - 1;
		ag->imps[i].a0 = &(ag->atoms[i1]);
		(ag->atoms[i1].nimps)++;

		atom2 = json_object_get(improper, "atom2");
		if (!json_is_integer(atom2)) {
			fprintf(stderr, "json atom2 is not integer for improper %zd in json_file %s\n", i, json_file);
		}
		int i2 = json_integer_value(atom2) - 1;
		ag->imps[i].a1 = &(ag->atoms[i2]);
		(ag->atoms[i2].nimps)++;

		atom3 = json_object_get(improper, "atom3");
		if (!json_is_integer(atom3)) {
			fprintf(stderr, "json atom3 is not integer for improper %zd in json_file %s\n", i, json_file);
		}
		int i3 = json_integer_value(atom3) - 1;
		ag->imps[i].a2 = &(ag->atoms[i3]);
		(ag->atoms[i3].nimps)++;

		atom4 = json_object_get(improper, "atom4");
		if (!json_is_integer(atom4)) {
			fprintf(stderr, "json atom4 is not integer for improper %zd in json_file %s\n", i, json_file);
		}
		int i4 = json_integer_value(atom4) - 1;
		ag->imps[i].a3 = &(ag->atoms[i4]);
		(ag->atoms[i4].nimps)++;

		phi = json_object_get(improper, "phi");
		if (!json_is_real(phi)) {
			fprintf(stderr, "json phi is not floating point for improper %zd in json_file %s\n", i, json_file);
		}
		ag->imps[i].psi0 = json_real_value(phi);

		spring_constant = json_object_get(improper, "spring_constant");
		if (!json_is_real(spring_constant)) {
			fprintf(stderr, "json spring_constant is not floating point for improper %zd in json_file %s\n", i, json_file);
		}
		ag->imps[i].k = json_real_value(spring_constant);
	}



	json_decref(base);

//allocate atom arrays of pointers to parameters
	for (size_t i = 0; i < natoms; i++) {
		int i1 = ag->atoms[i].nbonds;
		ag->atoms[i].bonds =
		    _mol_calloc(i1, sizeof(struct atombond *));
		ag->atoms[i].nbonds = 0;
		i1 = ag->atoms[i].nangs;
		ag->atoms[i].angs =
		    _mol_calloc(i1, sizeof(struct atomangle *));
		ag->atoms[i].nangs = 0;
		i1 = ag->atoms[i].ntors;
		ag->atoms[i].tors =
		    _mol_calloc(i1, sizeof(struct atomtorsion *));
		ag->atoms[i].ntors = 0;
		i1 = ag->atoms[i].nimps;
		ag->atoms[i].imps =
		    _mol_calloc(i1, sizeof(struct atomimproper *));
		ag->atoms[i].nimps = 0;
	}
	struct atom *atm;
//fill bonds
	for (int i = 0; i < ag->nbonds; i++) {
		atm = ag->bonds[i].a0;
		atm->bonds[(atm->nbonds)++] = &(ag->bonds[i]);
		atm = ag->bonds[i].a1;
		atm->bonds[(atm->nbonds)++] = &(ag->bonds[i]);
	}
//fill angles
	for (int i = 0; i < ag->nangs; i++) {
		atm = ag->angs[i].a0;
		atm->angs[(atm->nangs)++] = &(ag->angs[i]);
		atm = ag->angs[i].a1;
		atm->angs[(atm->nangs)++] = &(ag->angs[i]);
		atm = ag->angs[i].a2;
		atm->angs[(atm->nangs)++] = &(ag->angs[i]);
	}
//fill torsions
	for (int i = 0; i < ag->ntors; i++) {
		atm = ag->tors[i].a0;
		atm->tors[(atm->ntors)++] = &(ag->tors[i]);
		atm = ag->tors[i].a1;
		atm->tors[(atm->ntors)++] = &(ag->tors[i]);
		atm = ag->tors[i].a2;
		atm->tors[(atm->ntors)++] = &(ag->tors[i]);
		atm = ag->tors[i].a3;
		atm->tors[(atm->ntors)++] = &(ag->tors[i]);
	}
//fill impropers
	for (int i = 0; i < ag->nimps; i++) {
		atm = ag->imps[i].a0;
		atm->imps[(atm->nimps)++] = &(ag->imps[i]);
		atm = ag->imps[i].a1;
		atm->imps[(atm->nimps)++] = &(ag->imps[i]);
		atm = ag->imps[i].a2;
		atm->imps[(atm->nimps)++] = &(ag->imps[i]);
		atm = ag->imps[i].a3;
		atm->imps[(atm->nimps)++] = &(ag->imps[i]);
	}
//atom indices in the group
	fill_ingrp(ag);

	ag->is_psf_read = true;

	// copy vals from deprecated to new data structures
	int atomsi;
	for (atomsi = 0; atomsi < ag->natoms; atomsi++) {
		_mol_atom_create_bond_indices(&ag->atoms[atomsi],
					      ag->atoms[atomsi].nbonds);
	}
	_mol_atom_group_copy_from_deprecated(ag);
}
#endif /* _NO_JANSSON_ */
