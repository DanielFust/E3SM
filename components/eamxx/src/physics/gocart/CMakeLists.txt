# Source files for GOCART process
set(GOCART_SRCS
  eamxx_gocart_process_interface.cpp
)

# Create library for GOCART
set(GOCART_LIBS "gocart")
add_library(gocart ${GOCART_SRCS})
target_compile_definitions(gocart PUBLIC EAMXX_HAS_GOCART)
target_link_libraries(gocart physics_share scream_share)
target_compile_options(gocart PUBLIC)

include(ScreamUtils)

# Add this to eamxx_physics
target_link_libraries(eamxx_physics INTERFACE gocart)
