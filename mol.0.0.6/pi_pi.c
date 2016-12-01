#include <jansson.h>
#include <string.h>
#include _MOL_INCLUDE_

#define PI_PI_DIST_CUTOFF 5.0
#define PI_PI_OFFSET_DIST_CUTOFF 3.0

//come up with our own ags, which generates cubes for pseudoatoms

//we can actually reuse many functions from nbenergy.c as they are not tightly bound

void add_residue_pseudoatoms(struct pi_pi_setup * pps, const struct atomgrp *ag)
{
	int current_pseudoatom = pps->npseudoatoms;
	for (int i = 0; i < ag->nres; i++) {
		int atomi = ag->iares[i];
		char * residue_name = ag->atoms[atomi].residue_name;
		if (strcmp(residue_name, "HIE") == 0
				|| strcmp(residue_name, "HIP") == 0
				|| strcmp(residue_name, "HIS") == 0
				|| strcmp(residue_name, "PHE") == 0
				|| strcmp(residue_name, "TYR") == 0) {
			pps->npseudoatoms += 1;
		}
		else if (strcmp(residue_name, "TRP") == 0) {
			pps->npseudoatoms += 2;
		}
	}

	//if no new rings found
	if (pps->npseudoatoms == current_pseudoatom) {
		return;
	}

	pps->pseudoatoms = realloc(pps->pseudoatoms, pps->npseudoatoms * sizeof(struct pi_pi_pseudoatom));

	for (int i = 0; i < ag->nres; i++) {
		int atomi = ag->iares[i];
		char * residue_name = ag->atoms[atomi].residue_name;
		if (strcmp(residue_name, "HIE") == 0) {
			struct pi_pi_pseudoatom * pa = &(pps->pseudoatoms[current_pseudoatom]);
			pa->natoms= 5;
			pa->atom_indices = malloc(5 * sizeof(int));
			for (int j = 0; j < 5; j++) {
				pa->atom_indices[j] = atomi+4+j;
			}
			current_pseudoatom += 1;
		}
		else if (strcmp(residue_name, "HIP") == 0) {
			struct pi_pi_pseudoatom * pa = &(pps->pseudoatoms[current_pseudoatom]);
			pa->natoms= 5;
			pa->atom_indices = malloc(5 * sizeof(int));
			for (int j = 0; j < 3; j++) {
				pa->atom_indices[j] = atomi+4+j;
			}
			for (int j = 3; j < 5; j++) {
				pa->atom_indices[j] = atomi+5+j;
			}
			current_pseudoatom += 1;
		}
		else if (strcmp(residue_name, "HIS") == 0) {
			struct pi_pi_pseudoatom * pa = &(pps->pseudoatoms[current_pseudoatom]);
			pa->natoms= 5;
			pa->atom_indices = malloc(5 * sizeof(int));
			for (int j = 0; j < 2; j++) {
				pa->atom_indices[j] = atomi+4+j;
			}
			for (int j = 2; j < 5; j++) {
				pa->atom_indices[j] = atomi+5+j;
			}
			current_pseudoatom += 1;
		}
		else if (strcmp(residue_name, "PHE") == 0) {
			struct pi_pi_pseudoatom * pa = &(pps->pseudoatoms[current_pseudoatom]);
			pa->natoms= 6;
			pa->atom_indices = malloc(6 * sizeof(int));
			for (int j = 0; j < 6; j++) {
				pa->atom_indices[j] = atomi+4+j;
			}
			current_pseudoatom += 1;
		}
		else if (strcmp(residue_name, "TYR") == 0) {
			struct pi_pi_pseudoatom * pa = &(pps->pseudoatoms[current_pseudoatom]);
			pa->natoms= 6;
			pa->atom_indices = malloc(6 * sizeof(int));
			for (int j = 0; j < 6; j++) {
				pa->atom_indices[j] = atomi+4+j;
			}
			current_pseudoatom += 1;
		}
		else if (strcmp(residue_name, "TRP") == 0) {
			struct pi_pi_pseudoatom * pa = &(pps->pseudoatoms[current_pseudoatom]);
			pa->natoms= 5;
			pa->atom_indices = malloc(5 * sizeof(int));
			for (int j = 0; j < 3; j++) {
				pa->atom_indices[j] = atomi+4+j;
			}
			for (int j = 3; j < 5; j++) {
				pa->atom_indices[j] = atomi+5+j;
			}
			current_pseudoatom += 1;

			pa = &(pps->pseudoatoms[current_pseudoatom]);
			pa->natoms= 6;
			pa->atom_indices = malloc(6 * sizeof(int));
			for (int j = 0; j < 3; j++) {
				pa->atom_indices[j] = atomi+5+j;
			}
			for (int j = 3; j < 6; j++) {
				pa->atom_indices[j] = atomi+8+j;
			}
			current_pseudoatom += 1;
		}
	}
}

void add_probe_pseudoatoms(struct pi_pi_setup * pps, const int atom_offset, const char *json_file)
{

	int current_pseudoatom = pps->npseudoatoms;

	json_error_t json_file_error;
	json_t *base = json_load_file(json_file, 0, &json_file_error);
	json_t *pseudoatoms;
	pseudoatoms = json_object_get(base, "pi_pi_pseudoatoms");

	if (pseudoatoms == NULL) {
		return;
	}

	if (!json_is_array(pseudoatoms)) {
		fprintf(stderr, "json pseudoatoms are not an array %s\n", json_file);
	}
	size_t npatoms = json_array_size(pseudoatoms);
	pps->npseudoatoms += npatoms;

	pps->pseudoatoms = realloc(pps->pseudoatoms, pps->npseudoatoms * sizeof(struct pi_pi_pseudoatom));

	for (size_t i = 0; i < npatoms; i++) {
		json_t *pseudoatom = json_array_get(pseudoatoms, i);
		if (!json_is_array(pseudoatom)) {
			fprintf(stderr, "json pseudoatom %zd is not an array %s\n", i, json_file);
		}

		size_t natoms = json_array_size(pseudoatom);
		struct pi_pi_pseudoatom * pa = &(pps->pseudoatoms[current_pseudoatom]);
		pa->natoms = natoms;
		pa->atom_indices = malloc(natoms * sizeof(int));

		json_t *atom_index;
		size_t j;
		json_array_foreach(pseudoatom, j, atom_index) {
			if (!json_is_integer(atom_index)) {
				fprintf(stderr, "json pseudoatom %zd atom %zd is not an integer %s\n", i, j, json_file);
			}
			pa->atom_indices[j] = atom_offset + json_integer_value(atom_index)-1;
		}
	}
}

static void update_centers(struct pi_pi_setup * pps, const struct atomgrp *ag)
{
	for (int i = 0; i < pps->npseudoatoms; i++) {
		struct pi_pi_pseudoatom *pa = &(pps->pseudoatoms[i]);
		pa->center.X = pa-> center.Y = pa->center.Z = 0.0;
		for (int j = 0; j < pa->natoms; j++) {
			int index = pa->atom_indices[j];
			pa->center.X += ag->atoms[index].X;
			pa->center.Y += ag->atoms[index].Y;
			pa->center.Z += ag->atoms[index].Z;
		}

		pa->center.X /= pa->natoms;
		pa->center.Y /= pa->natoms;
		pa->center.Z /= pa->natoms;
	}
}

static double pi_pi_offset_distance(struct pi_pi_pseudoatom *par,
		struct pi_pi_pseudoatom *pal, struct mol_vector3 *plane_r)
{
	double rec_d = (plane_r->X * par->center.X) + (plane_r->Y * par->center.Y) +
			(plane_r->Z * par->center.Z);

	double lig_t = rec_d - (plane_r->X * pal->center.X) - (plane_r->Y * pal->center.Y) -  (plane_r->Z * pal->center.Z);
	
	double offset1_x = par->center.X - (pal->center.X+(lig_t*plane_r->X));
	double offset1_y = par->center.Y - (pal->center.Y+(lig_t*plane_r->Y));
	double offset1_z = par->center.Z - (pal->center.Z+(lig_t*plane_r->Z));
	    
	double offset = sqrt(_mol_sq(offset1_x) + _mol_sq(offset1_y) + _mol_sq(offset1_z))
;
	return offset;
}

static double pi_pi_pairwise_eng(struct atomgrp *ag, struct pi_pi_pseudoatom *par,
		struct pi_pi_pseudoatom *pal, double weight)
{
	double energy;

	//receptor
	struct mol_vector3 ar[4];
	struct mol_vector3 plane_r;

	//ligand
	struct mol_vector3 al[4];
	struct mol_vector3 plane_l;

	for (int i =0; i < 4; i++) {
		int atomr = par->atom_indices[i];
		int atoml = pal->atom_indices[i];

		ar[i].X = ag->atoms[atomr].X;
		ar[i].Y = ag->atoms[atomr].Y;
		ar[i].Z = ag->atoms[atomr].Z;

		al[i].X = ag->atoms[atoml].X;
		al[i].Y = ag->atoms[atoml].Y;
		al[i].Z = ag->atoms[atoml].Z;
	}
	//ar[0].X += ar[1].X;
	//ar[0].Y += ar[1].Y;
	//ar[0].Z += ar[1].Z;
	//ar[1].X += ar[2].X;
	//ar[1].Y += ar[2].Y;
	//ar[1].Z += ar[2].Z;
	//ar[2].X += ar[3].X;
	//ar[2].Y += ar[3].Y;
	//ar[2].Z += ar[3].Z;

	//al[0].X += al[1].X;
	//al[0].Y += al[1].Y;
	//al[0].Z += al[1].Z;
	//al[1].X += al[2].X;
	//al[1].Y += al[2].Y;
	//al[1].Z += al[2].Z;
	//al[2].X += al[3].X;
	//al[2].Y += al[3].Y;
	//al[2].Z += al[3].Z;

	//ar[0].X -= ar[1].X;
	//ar[0].Y -= ar[1].Y;
	//ar[0].Z -= ar[1].Z;
	//ar[1].X -= ar[2].X;
	//ar[1].Y -= ar[2].Y;
	//ar[1].Z -= ar[2].Z;

	//al[0].X -= al[1].X;
	//al[0].Y -= al[1].Y;
	//al[0].Z -= al[1].Z;
	//al[1].X -= al[2].X;
	//al[1].Y -= al[2].Y;
	//al[1].Z -= al[2].Z;
	ar[0].X -= ar[2].X;
	ar[0].Y -= ar[2].Y;
	ar[0].Z -= ar[2].Z;
	ar[1].X -= ar[3].X;
	ar[1].Y -= ar[3].Y;
	ar[1].Z -= ar[3].Z;

	al[0].X -= al[2].X;
	al[0].Y -= al[2].Y;
	al[0].Z -= al[2].Z;
	al[1].X -= al[3].X;
	al[1].Y -= al[3].Y;
	al[1].Z -= al[3].Z;

	//(yn,n+1*zn+1,n+2 - zn,n+1*yn+1,n+2; zn,n+1*xn+1,n+2 - xn,n+1*zn+1,n+2; xn,n+1*yn+1,n+2 - yn,n+1*xn+1,n+2)
	plane_r.X = ar[0].Y*ar[1].Z - ar[0].Z*ar[1].Y;
	plane_r.Y = ar[0].Z*ar[1].X - ar[0].X*ar[1].Z;
	plane_r.Z = ar[0].X*ar[1].Y - ar[0].Y*ar[1].X;

	plane_l.X = al[0].Y*al[1].Z - al[0].Z*al[1].Y;
	plane_l.Y = al[0].Z*al[1].X - al[0].X*al[1].Z;
	plane_l.Z = al[0].X*al[1].Y - al[0].Y*al[1].X;

	double rlen = sqrt(_mol_sq(plane_r.X) + _mol_sq(plane_r.Y) + _mol_sq(plane_r.Z));
	plane_r.X /= rlen;
	plane_r.Y /= rlen;
	plane_r.Z /= rlen;

	double offset = pi_pi_offset_distance(par, pal, &plane_r);
	//fprintf(stderr, "Offset: %f\n", offset);

	if (offset > PI_PI_OFFSET_DIST_CUTOFF) {
		return 0.0;
	}

//	E = f(x)/h(x)¬
//	f(x) = abs(g(x))¬
//	E = abs(g(x))/h(x)¬

	double llen = sqrt(_mol_sq(plane_l.X) + _mol_sq(plane_l.Y) + _mol_sq(plane_l.Z));
	double h_x = llen;

	double g_x = plane_r.X*plane_l.X + plane_r.Y*plane_l.Y + plane_r.Z * plane_l.Z;
	double f_x = fabs(g_x);
	energy = f_x / llen;

//	f'(x) = g(x)*g'(x)/abs(g(x))¬
//	f'(x) = g(x)*g'(x)/f(x)¬
//	E' = (f'(x)h(x) - f(x)h'(x)) / (h(x)**2)¬
//	¬
//	dg/dxn = -zn+1,n+2 * Ys + yn+1,n+2 * Zs¬
//	dg/dyn = -xn+1,n+2 * Zs + zn+1,n+2 * Xs¬
//	dg/dzn = -yn+1,n+2 * Xs + xn+1,n+2 * Ys¬
//	¬
//	dh/dxn = ((yn+1,n+2 - zn+1,n+2) * xn) / h¬
//	dh/dyn = ((zn+1,n+2 - xn+1,n+2) * yn) / h¬
//	dh/dzn = ((xn+1,n+2 - yn+1,n+2) * zn) / h¬

	for (int i = 0; i < pal->natoms; i++) {
		int i1 = i + 1;
		if (i1 == pal->natoms) {
			i1 = 0;
		}
		int i2 = i1 + 1;
		if (i2 == pal->natoms) {
			i2 = 0;
		}
		int i3 = i2 + 1;
		if (i3 == pal->natoms) {
			i3 = 0;
		}
		int n = pal->atom_indices[i];
		int n1 = pal->atom_indices[i1];
		int n2 = pal->atom_indices[i2];
		int n3 = pal->atom_indices[i3];

		struct atom *an = (&ag->atoms[n]);
		struct atom *an1 = (&ag->atoms[n1]);
		struct atom *an2 = (&ag->atoms[n2]);
		struct atom *an3 = (&ag->atoms[n3]);

		double Xn = an->X;
		double Yn = an->Y;
		double Zn = an->Z;

		double Xn1n2n3 = an1->X - an3->X;
		double Yn1n2n3 = an1->Y - an3->Y;
		double Zn1n2n3 = an1->Z - an3->Z;

		double Xn0n1n2 = an->X - an2->X;
		double Yn0n1n2 = an->Y - an2->Y;
		double Zn0n1n2 = an->Z - an2->Z;

		struct mol_vector3 df;
		struct mol_vector3 dg;
		struct mol_vector3 dh;
		struct mol_vector3 dE;

		dg.X = (-Zn1n2n3 * plane_r.Y) + (Yn1n2n3 * plane_r.Z);
		dg.Y = (-Xn1n2n3 * plane_r.Z) + (Zn1n2n3 * plane_r.X);
		dg.Z = (-Yn1n2n3 * plane_r.X) + (Xn1n2n3 * plane_r.Y);

		dh.X = (_mol_sq(Zn1n2n3) + _mol_sq(Yn1n2n3)) * Xn;
		dh.X += -Zn1n2n3 * (Zn0n1n2 * Xn1n2n3 + an2->X * Zn1n2n3);
		dh.X += Yn1n2n3 * (-an2->X * Yn1n2n3 - Yn0n1n2 * Xn1n2n3);
		dh.X /= h_x;

		dh.Y = (_mol_sq(Xn1n2n3) + _mol_sq(Zn1n2n3)) * Yn;
		dh.Y += -Xn1n2n3 * (Xn0n1n2 * Yn1n2n3 + an2->Y * Xn1n2n3);
		dh.Y += Zn1n2n3 * (-an2->Y * Zn1n2n3 - Zn0n1n2 * Yn1n2n3);
		dh.Y /= h_x;

		dh.Z = (_mol_sq(Yn1n2n3) + _mol_sq(Xn1n2n3)) * Zn;
		dh.Z += -Yn1n2n3 * (Yn0n1n2 * Zn1n2n3 + an2->Z * Yn1n2n3);
		dh.Z += Xn1n2n3 * (-an2->Z * Xn1n2n3 - Xn0n1n2 * Zn1n2n3);
		dh.Z /= h_x;

		df.X = g_x*dg.X/f_x;
		df.Y = g_x*dg.Y/f_x;
		df.Z = g_x*dg.Z/f_x;

		dE.X = ((df.X * h_x) - (f_x * dh.X))/(h_x*h_x);
		dE.Y = ((df.Y * h_x) - (f_x * dh.Y))/(h_x*h_x);
		dE.Z = ((df.Z * h_x) - (f_x * dh.Z))/(h_x*h_x);

		an->GX += dE.X * weight;
		an->GY += dE.Y * weight;
		an->GZ += dE.Z * weight;
	}

	//fprintf(stderr, "%f\n", energy);
	return -energy*weight;
}

void pi_pi_eng(struct atomgrp *ag, double *energy, struct pi_pi_setup *pps_rec,
		struct pi_pi_setup *pps_lig, double weight)
{
	update_centers(pps_rec, ag);
	update_centers(pps_lig, ag);

	for (int i = 0; i < pps_rec->npseudoatoms; i++) {
		struct pi_pi_pseudoatom *par = &(pps_rec->pseudoatoms[i]);
		for (int j = 0; j < pps_lig->npseudoatoms; j++) {
			struct pi_pi_pseudoatom *pal = &(pps_lig->pseudoatoms[j]);
			double dsq = _mol_sq(par->center.X - pal->center.X) +
				_mol_sq(par->center.Y - pal->center.Y) +
				_mol_sq(par->center.Z - pal->center.Z);

			if (dsq < (PI_PI_DIST_CUTOFF * PI_PI_DIST_CUTOFF)) {
				//fprintf(stderr, "%f\n", dsq);
				(*energy) += pi_pi_pairwise_eng(ag, par, pal, weight);
			}
		}
	}

}
