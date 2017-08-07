# Setup optional dependencies and features
# alters these things
#
#       MOLSTURM_DEPENDENCIES			everyone needs these libraries
#       MOLSTURM_DEPENDENCIES_DEBUG		debug mode needs these extras
#       MOLSTURM_DEPENDENCIES_RELEASE		release mode needs these extras
#       MOLSTURM_DEPENDENCIES_TEST		tests need these extra libraries
#

####################
#-- C++ standard --#
####################
if (DRB_HAS_CXX14_SUPPORT)
	message(STATUS "Detected C++14 support: Setting MOLSTURM_HAVE_CXX14")
	set(MOLSTURM_HAVE_CXX14 ON)
endif()
if (DRB_HAS_CXX17_SUPPORT)
	message(STATUS "Detected C++17 support: Setting MOLSTURM_HAVE_CXX17")
	set(MOLSTURM_HAVE_CXX17 ON)
endif()
