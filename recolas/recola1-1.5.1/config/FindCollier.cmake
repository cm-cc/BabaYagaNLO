# File: FindCollier.cmake
# Author: Jean-Nicolas Lang
# Description: Find collier library using different approaches
# Last Modified: August 04, 2017

function(find_collier CLL_LIB_PATH CLL_INC_PATH)
  # Option 1: collier_path given
  if (collier_path)
    set (CLL_LIB_PATH ${collier_path})
    set (CLL_INC_PATH ${collier_path}/modules)
    unset(collier_path CACHE)

    find_library(COLLIER_LIBRARY
                 NAMES libcollier collier
                 HINTS "${CLL_LIB_PATH}"
    )
    if (COLLIER_LIBRARY)
      find_path(COLLIER_INCLUDE
                NAMES collier.mod 
                HINTS "${CLL_LIB_PATH}/modules" "${CLL_LIB_PATH}/include"
      )
      if (COLLIER_INCLUDE)
        set (CLL_LIB_PATH ${COLLIER_LIBRARY} PARENT_SCOPE)
        set (CLL_INC_PATH ${COLLIER_INCLUDE} PARENT_SCOPE)
        unset(COLLIER_LIBRARY CACHE)
        unset(COLLIER_INCLUDE CACHE)
      else()
        message(WARNING "Found collier library in ${CLL_LIB_PATH}, but no header files.")
        unset(COLLIER_LIBRARY CACHE)
        message("COLLIER_LIBRARY ${COLLIER_LIBRARY}")
      endif()
    endif()
  # Option 2: COLLIER_LIBRARY_DIR & COLLIER_INCLUDE_DIR given
  elseif (DEFINED ENV{COLLIER_INCLUDE_DIR} AND DEFINED ENV{COLLIER_LIBRARY_DIR})
    set (CLL_LIB_PATH $ENV{COLLIER_LIBRARY_DIR})
    set (CLL_INC_PATH $ENV{COLLIER_INCLUDE_DIR})
    find_library(COLLIER_LIBRARY
                 NAMES libcollier collier
                 HINTS "${CLL_LIB_PATH}"
    )
    if (COLLIER_LIBRARY)
      find_path(COLLIER_INCLUDE
                NAMES collier.mod 
                HINTS "${CLL_INC_PATH}"
      )
      if (COLLIER_INCLUDE)
        set (CLL_LIB_PATH ${COLLIER_LIBRARY} PARENT_SCOPE)
        set (CLL_INC_PATH ${COLLIER_INCLUDE} PARENT_SCOPE)
        unset(COLLIER_LIBRARY CACHE)
        unset(COLLIER_INCLUDE CACHE)
      else()
        unset(COLLIER_LIBRARY)
      endif()
    endif()
  # Option 3: Searching for a collier directory RELATIVE TO the model file path
  else()
    get_filename_component(CLL_LIB_PATH ${PROJECT_SOURCE_DIR} DIRECTORY)
    file(GLOB POTENTIAL_COLLIER_PATHS "${CLL_LIB_PATH}/COLLIER-*")
    foreach (CLL_LIB_PATH ${POTENTIAL_COLLIER_PATHS})
      if(IS_DIRECTORY "${CLL_LIB_PATH}")
        message("Searching collier library in path: ${CLL_LIB_PATH}")
        find_library(COLLIER_LIBRARY
                     NAMES libcollier collier
                     HINTS "${CLL_LIB_PATH}"
        )
        if (COLLIER_LIBRARY)
          find_path(COLLIER_INCLUDE
                    NAMES collier.mod 
                    HINTS "${CLL_LIB_PATH}/modules" "${CLL_LIB_PATH}/include"
          )
          if (COLLIER_INCLUDE)
            set (CLL_LIB_PATH ${COLLIER_LIBRARY} PARENT_SCOPE)
            set (CLL_INC_PATH ${COLLIER_INCLUDE} PARENT_SCOPE)
            unset(COLLIER_LIBRARY CACHE)
            unset(COLLIER_INCLUDE CACHE)
            break()
          else()
            message(WARNING "Found collier library in ${CLL_LIB_PATH}, but no header files.")
            unset(COLLIER_LIBRARY CACHE)
            message("COLLIER_LIBRARY ${COLLIER_LIBRARY}")
          endif()
        endif()
      endif()
    endforeach()
  endif()
endfunction()
