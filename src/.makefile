#------------------------------------------------------------------------------
# AFC -- MAKEFILE
#------------------------------------------------------------------------------
#
# Tobias Holzmann
# March 2020
#
# Description
#     This makefile includes all messages that are printed during compilation
#
#------------------------------------------------------------------------------

.PHONEY = buildThermoKinetics pre linkHeaders libThermoKinetics

PRE_MSG=\
	$(info ) \
	$(info ************************* AFC MAKE ****************************) \
	$(info ) \
	$(info c-o Starting compiling the libraries) \
	$(info c-o Build by Tobias Holzmann) \
	$(info c-o Makefile update 02.03.2020) \
	$(info ) \
	$(info ***************************************************************) \
	$(info )



POST_MSG=\
	$(info ) \
	$(info ) \
	$(info ***************************************************************) \
	$(info ) \
	$(info c-o Compiling of the libraries done) \
	$(info c-o You are ready to use the flamelet creator or the libraries) \
	$(info c-o If you have any problems comment on git or write an email) \
	$(info c-o Makefile update 02.03.2020) \
	$(info ) \
	$(info ***************************************************************)


CLEAN_MSG=\
	$(info ) \
	$(info ************************* AFC MAKE ****************************) \
	$(info ) \
	$(info c-o Cleaning the libraries and folder) \
	$(info ) \
	$(info ***************************************************************) \
	$(info )

#------------------------------------------------------------------------------
