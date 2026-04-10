# Common flags shared across platform configs.
# Included by platform-specific .mk files if needed.

# Debug build overrides
ifdef DEBUG
  ifeq ($(OS),Windows_NT)
    CXXFLAGS += /Od /Zi /DDEBUG
  else
    CXXFLAGS += -O0 -g -DDEBUG
  endif
endif
