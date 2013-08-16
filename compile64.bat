del mol.0.0.6\*.obj mol64.lib
cl /Ox /TC /MP4 /nologo /D:_NO_JANSSON_ /c /Fo.\mol.0.0.6\ .\mol.0.0.6\atom.c .\mol.0.0.6\atom_group.c .\mol.0.0.6\_atom_group_copy_from_deprecated.c .\mol.0.0.6\benergy.c .\mol.0.0.6\bond.c .\mol.0.0.6\compare.c .\mol.0.0.6\energy.c .\mol.0.0.6\gbsa.c .\mol.0.0.6\hbond.c .\mol.0.0.6\icharmm.c .\mol.0.0.6\io.c .\mol.0.0.6\lbfgs.c .\mol.0.0.6\matrix.c .\mol.0.0.6\mem.c .\mol.0.0.6\mol2.c .\mol.0.0.6\ms.c .\mol.0.0.6\myhelpers.c .\mol.0.0.6\nbenergy.c .\mol.0.0.6\octree.c .\mol.0.0.6\pdb.c .\mol.0.0.6\potential.c .\mol.0.0.6\prms.c .\mol.0.0.6\protein.c .\mol.0.0.6\rmsd.c .\mol.0.0.6\rotamer.c .\mol.0.0.6\sasa.c .\mol.0.0.6\sdf.c .\mol.0.0.6\subag.c .\mol.0.0.6\version.c

lib /NOLOGO /OUT:mol64.lib mol.0.0.6\*.obj
