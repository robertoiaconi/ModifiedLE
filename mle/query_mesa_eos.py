import pyMesaUtils as pym
import numpy as np
import os

#check if the mesa directory is set
if "MESA_DIR" not in os.environ:
    raise ValueError("Must set MESA_DIR environment variable or install MESA")

#initialise necessary MESA modules
eos_lib, eos_def = pym.loadMod("eos")
const_lib, const_def = pym.loadMod("const")
if pym.MESA_VERSION < 12608:
	crlibm_lib, _ = pym.loadMod("crlibm")
	crlibm_lib.crlibm_init()
else:
	crlibm_lib, _ = pym.loadMod("math")
	crlibm_lib.math_init()
chem_lib, chem_def = pym.loadMod("chem")
net_lib, net_def = pym.loadMod("net")
rates_lib, rates_def = pym.loadMod("rates")

ierr=0	
const_lib.const_init(pym.MESA_DIR,ierr)
chem_lib.chem_init('isotopes.data',ierr) #chem lib needs const and crlibm libs to work appropriately

ierr=0
if pym.MESA_VERSION >= 10000:
     #Function sig changed
     rates_lib.rates_init('reactions.list','jina_reaclib_results_20130213default2','rate_tables',False,False,'','','',ierr)
else:
     rates_lib.rates_init('reactions.list','jina_reaclib_results_20130213default2','rate_tables',False,'','','',ierr)
net_lib.net_init(ierr) #net lib needs the rates lib to work appropriately
eos_lib.eos_init('mesa','','','',True,ierr)

ierr=0              
eos_handle = eos_lib.alloc_eos_handle(ierr)

#queries the MESA EoS for the gas internal energy, from density and pressure and solving for temperature with an iteration method
def query_mesa_eos_solveT(Z,X,abar,zbar,species,chem_id,net_iso,mfrac,log10Rho,log10P,logT_guess):
	#gets elements chemistries and networks to include when calling the EoS
	chemistry = np.array([])
	for i,value in enumerate(chem_id):
		if value not in dir(chem_def):
    			raise ValueError(value+" not in the list of valid MESA elements.")

		chemistry = np.append(chemistry, eval('chem_def.'+value+'.get()'))

	reactions = np.array([])
	for i,value in enumerate(net_iso):
		if value not in dir(net_lib):
    			raise ValueError(value+" not in the list of valid MESA net isotopes.")

		reactions = np.append(reactions, eval('net_lib.'+value+'.get()'))
	
	#query mesa eos, see a detailed explanation in the definition of the function "eosDT_get_T_given_Ptotal" in $MESA_DIR/eos/public/eos_lib.f90
	 #in (some quantities hard coded here according to thair values in MESA)
	#handle, Z, X, abar, zbar, species, chem_id, net_iso, xa, logRho, logP, logP_tol, logT_tol, max_iter, logT_guess, logT_bnd1, logT_bnd2, logP_at_bnd1, logP_at_bnd2
	logT_tol = 1.e-8
	logP_tol = 1.e-8
	max_iter = 100
	logT_bnd1 = 9e99
	logT_bnd2 = 9e99
	logP_at_bnd1 = 9e99
	logP_at_bnd2 = 9e99

	 #out
	logT_result = 0.0
	res = np.zeros(eos_def.num_eos_basic_results.get())
	d_dlnRho_const_T = np.zeros(eos_def.num_eos_basic_results.get())
	d_dlnT_const_Rho = np.zeros(eos_def.num_eos_basic_results.get())
	d_dabar_const_TRho = np.zeros(eos_def.num_eos_basic_results.get())
	d_dzbar_const_TRho = np.zeros(eos_def.num_eos_basic_results.get())
	eos_calls = 0
	ierr = 0
	
     #query
	eos_res = eos_lib.eosdt_get_t_given_ptotal(eos_handle, Z, X, abar, zbar, species, chemistry, reactions, mfrac, log10Rho, log10P,\
                                                   logT_tol, logP_tol, max_iter,\
                                                   logT_guess,\
                                                   logT_bnd1, logT_bnd2, logP_at_bnd1, logP_at_bnd2,\
                                                   logT_result, res, d_dlnRho_const_T, d_dlnT_const_Rho, d_dabar_const_TRho, d_dzbar_const_TRho, eos_calls, ierr)

	#save necessary quantities
	res = eos_res["res"]
	res_temp = eos_res["logt_result"] #dictionary keys are case sensitive	

	i_lnE = eos_def.i_lnE.get() - 1
	IE = np.exp(res[i_lnE])
#	print("Internal Energy from MESA EoS: ", IE, " erg/g")

	i_Pgas = eos_def.i_lnpgas.get() - 1
	PGAS = np.exp(res[i_Pgas])
#	print("Gas pressure from MESA EoS: ", PGAS, " erg/cm^3")

	T = 10.**res_temp
#	print("Gas temperature from MESA EoS: ", PGAS, " erg/cm^3")

	return IE, PGAS, T

