include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCE_FILES    wdenmath.c
                    ebayesthresh.c
                    wdenoise.c
                    )

set(HEADER_FILES    wdenmath.h
                    ebayesthresh.h
                    )

add_library(wdenoiselib STATIC ${SOURCE_FILES} ${HEADER_FILES})

set_property(TARGET wdenoiselib PROPERTY FOLDER "lib")

target_include_directories(wdenoiselib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}../header)
