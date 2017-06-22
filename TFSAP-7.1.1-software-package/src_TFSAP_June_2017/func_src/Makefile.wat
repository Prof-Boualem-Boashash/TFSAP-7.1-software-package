#############################################################################
#
# Makefile for TFSA (Win version)
#
# Tested with Open Watcom C++ compiler
# 
# Do not edit!
#
#############################################################################

# Inference rule for compiling C files into O files
.c.obj:
	$(CC) $(CFLAGS) $<  /fo..\$@
