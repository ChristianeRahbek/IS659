function(add_subproject SUBPROJECTS)
  foreach(sp ${SUBPROJECTS})
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/data  ${CMAKE_SOURCE_DIR}/${sp}/data
                    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/setup  ${CMAKE_SOURCE_DIR}/${sp}/setup)
    add_subdirectory(${CMAKE_SOURCE_DIR}/${sp})
  endforeach()
endfunction()