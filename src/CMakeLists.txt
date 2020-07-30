dune_symlink_to_source_files(FILES domain.msh plane-elasticity.input)
add_executable("plane-elasticity" plane-elasticity.cc driver.hh operator.hh bctype.hh plane-elasticity.input)
target_link_dune_default_libraries("plane-elasticity")

