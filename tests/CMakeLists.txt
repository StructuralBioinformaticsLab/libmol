include_directories(..)

include_directories(${CHECK_INCLUDE_DIRS})

# Add tests
add_executable(test_mol_atom test_mol_atom.c)
target_link_libraries(test_mol_atom
  ${CHECK_LIBRARIES}
  mol.${libmol_version})
add_executable(test_mol_pdb test_mol_pdb.c)
target_link_libraries(test_mol_pdb
  ${CHECK_LIBRARIES}
  mol.${libmol_version})
add_executable(test_benergy test_benergy.c)
target_link_libraries(test_benergy
  ${CHECK_LIBRARIES}
  mol.${libmol_version} m)

# Configure data files
file(GLOB test_files "${CMAKE_CURRENT_SOURCE_DIR}/data/*")
foreach(filepath ${test_files})
  get_filename_component(filename ${filepath} NAME)
  configure_file(${filepath} ${CMAKE_CURRENT_BINARY_DIR}/${filename} COPYONLY)
endforeach()

# Register tests
add_test(test_mol_atom ${CMAKE_CURRENT_BINARY_DIR}/test_mol_atom)
add_test(test_mol_pdb ${CMAKE_CURRENT_BINARY_DIR}/test_mol_pdb)
add_test(test_benergy ${CMAKE_CURRENT_BINARY_DIR}/test_benergy)
