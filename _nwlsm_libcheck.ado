*! version 0.2.0  19nov2026 nwlsm
* Copyright (C) 2026 Thomas U. Grund. Part of nwlsm; GNU GPL-3 (see LICENSE).
capture program drop _nwlsm_libcheck
program _nwlsm_libcheck
	version 14
	capture mata: mata which LsmFit()
	if _rc {
		capture mata: mata mlib index
		capture mata: mata which LsmFit()
		if _rc {
			di as error "{p 0 0 2}nwlsm: the Mata library lnwlsm.mlib was not found. Reinstall the package with {bf:net install nwlsm, replace}.{p_end}"
			exit 3499
		}
	}
end
