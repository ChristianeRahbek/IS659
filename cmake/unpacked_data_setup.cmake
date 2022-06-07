set(unpacked_data_dir "/mnt/rust/data/is659/unpacked")
if(NOT EXISTS ${unpacked_data_dir})
  message(FATAL_ERROR "The path ${unpacked_data_dir} does not exist. Please specify the correct path to the unpacked "
      "data of the experiment, i.e. pixie2ausa-converted .root files with names 'Run###.root', "
      "in the file ${CMAKE_CURRENT_LIST_FILE} close to where this error is thrown.")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/data
                COMMAND ${CMAKE_COMMAND} -E create_symlink ${unpacked_data_dir}  ${CMAKE_SOURCE_DIR}/data/unpacked
                COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/data  ${CMAKE_SOURCE_DIR}/analysis/data
                COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/data  ${CMAKE_SOURCE_DIR}/calibration/data)