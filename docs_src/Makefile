# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = ..\docs

.PHONY: html Makefile

html: Makefile
	@rm -r "$(BUILDDIR)"/*
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@rm -r "$(BUILDDIR)"/doctrees
	@mv "$(BUILDDIR)"/html/* "$(BUILDDIR)"/
	@rm -r "$(BUILDDIR)"/html
