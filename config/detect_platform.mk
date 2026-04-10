# Auto-detect OS, architecture, and available compilers.
# Sets PLATFORM variable and includes the appropriate platform config.

ifeq ($(OS),Windows_NT)
    PLATFORM := windows_x86_64
    include config/windows_x86_64.mk
else
    UNAME_S := $(shell uname -s)
    UNAME_M := $(shell uname -m)
    ifeq ($(UNAME_S),Linux)
        ifeq ($(UNAME_M),x86_64)
            PLATFORM := linux_x86_64
        else
            PLATFORM := linux_$(UNAME_M)
        endif
    else ifeq ($(UNAME_S),Darwin)
        PLATFORM := macos_$(UNAME_M)
    else
        PLATFORM := unknown
    endif
    include config/linux_x86_64.mk
endif
