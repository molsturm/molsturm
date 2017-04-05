set(SRC main.cc parse_args.cc)

# TODO Only compile release version of hf:
# add_executable(hf parse_args.cc ${SRC})
# setup_example_target(hf)

# XXX: For now compile hf with DRB
set(TARGET hf)
drb_add_executable(${TARGET} 
	DBGSUFFIX  "_DEBUG"
	FILES ${SRC}
)

drb_target_compile_definitions(ALL     ${TARGET}  PRIVATE ${MOLSTURM_DEFINITIONS})
drb_target_compile_definitions(DEBUG   ${TARGET}  PRIVATE ${MOLSTURM_DEFINITIONS_DEBUG})
drb_target_compile_definitions(RELEASE ${TARGET}  PRIVATE ${MOLSTURM_DEFINITIONS_RELEASE})

drb_target_link_libraries(ALL     ${TARGET}  ${MOLSTURM_DEPENDENCIES})
drb_target_link_libraries(DEBUG   ${TARGET}  ${molsturm_DEBUG_TARGET} ${MOLSTURM_DEPENDENCIES_DEBUG})
drb_target_link_libraries(RELEASE ${TARGET}  ${molsturm_RELEASE_TARGET} ${MOLSTURM_DEPENDENCIES_RELEASE})
