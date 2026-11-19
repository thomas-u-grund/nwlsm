*! version 0.5.0  19nov2026 nwlsm
* Copyright (C) 2026 Thomas U. Grund. Part of nwlsm; GNU GPL-3 (see LICENSE).
capture program drop nwlsm_version
program nwlsm_version, rclass
	version 14
	local v "0.1.0"
	local d "19nov2026"
	di as text "nwlsm version `v', `d'"
	di as text "  status                  : development (base Euclidean model)"
	di as text "  license                 : GNU GPL-3"
	return local version "`v'"
	return local date "`d'"
end
