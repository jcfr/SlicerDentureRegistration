cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(DentureRegistration)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/gongmingjun/SlicerDentureRegistration#readme")
set(EXTENSION_CATEGORY "Registration")
set(EXTENSION_CONTRIBUTORS "Xiaojun Chen (SJTU), Mingjun Gong (SJTU), Yueang Liu (SJTU)")
set(EXTENSION_DESCRIPTION "This is an extension module to realize the registration of patients' CT images and the corresponding denture models.")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/gongmingjun/SlicerDentureRegistration/master/icon.png")
set(EXTENSION_SCREENSHOTURLS "https://raw.githubusercontent.com/gongmingjun/SlicerDentureRegistration/master/DentureRegistration.jpg")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(Repair)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
