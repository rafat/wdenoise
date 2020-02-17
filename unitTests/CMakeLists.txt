
set(SOURCE_FILES
		unitTest.c
	)

set(HEADER_FILES
        unitTest.h
    )


add_executable(wdenoiseTests ${SOURCE_FILES} )

add_test(NAME wdenoiseTests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test  COMMAND wdenoiseTests)

add_dependencies(wdenoiseTests wdenoiselib wavelib)
target_link_libraries(wdenoiseTests wdenoiselib wavelib m)

target_include_directories(wdenoiseTests PUBLIC
							${CMAKE_SOURCE_DIR}/../src
                            )


install(TARGETS wdenoiseTests
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION tests
    ARCHIVE DESTINATION tests
    )
