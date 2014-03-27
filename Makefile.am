SUBDIRS = \
	src \
        test

# doxygen exists
if DOX_SET 
# make doxygen docs (as defined by Doxyfile)

# currently Doxyfile if set to make html, latex and man pages
# this could be changed, or separate config (with .cfg extension)
# file could be made for each
DOX_CFG_FILE=Doxyfile
DOC_FILES = \
  $(DOX_CFG_FILE) \
  src/ampmarginaliser.dox \
  src/ampmarginaliser.h \
  src/ampmarginaliser.c

doc: $(DOC_FILES)
	doxygen $(DOX_CFG_FILE)
endif

