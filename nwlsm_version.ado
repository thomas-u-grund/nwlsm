*! version 1.0.1  20nov2026 nwlsm
* Copyright (C) 2026 Thomas U. Grund. Part of nwlsm; GNU GPL-3 (see LICENSE).
capture program drop nwlsm_version
program nwlsm_version, rclass
	version 14
	local __pkgver "1.0.1"
	local __pkgdate "20nov2026"
	di as text "nwlsm version `__pkgver', `__pkgdate'"
	di as text "  status                  : latent space models (Euclidean and bilinear, valued, clusters, influence)"
	di as text "  license                 : GNU GPL-3"
	return local version "`__pkgver'"
	return local date "`__pkgdate'"
end
