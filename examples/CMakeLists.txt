add_executable(denoiseexamples denoiseexamples.c test2d.c)

target_link_libraries(denoiseexamples wdenoiselib wavelib)

add_executable(wdenoiseexample wdenoiseexample.c test2d.c)

target_link_libraries(wdenoiseexample wdenoiselib wavelib)

add_executable(ebayesdenoiseexample ebayesdenoiseexample.c test2d.c)

target_link_libraries(ebayesdenoiseexample wdenoiselib wavelib)

add_executable(imdenoiseexample imdenoiseexample.c test2d.c)

target_link_libraries(imdenoiseexample wdenoiselib wavelib)


if(UNIX)
	target_link_libraries(denoiseexamples m)
	target_link_libraries(wdenoiseexample m)
	target_link_libraries(ebayesdenoiseexample m)
  target_link_libraries(imdenoiseexample m)
endif()

set_target_properties(denoiseexamples wdenoiseexample ebayesdenoiseexample imdenoiseexample
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/examples"
        )
