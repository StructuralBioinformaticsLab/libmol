include_directories(..)

include_directories(${CHECK_INCLUDE_DIRS})

add_executable(test_mol_atom test_mol_atom.c)
target_link_libraries(test_mol_atom
  ${CHECK_LIBRARIES}
  mol.${libmol_version})
add_test(test_mol_atom ${CMAKE_CURRENT_BINARY_DIR}/test_mol_atom)
