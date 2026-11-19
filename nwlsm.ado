*! version 0.1.0  19nov2026 nwlsm
* Copyright (C) 2026 Thomas U. Grund
* This file is part of nwlsm, a Stata package for latent space models for networks. It follows the latent space model of
* Hoff, Raftery and Handcock (2002) as implemented in the R package latentnet (Krivitsky and Handcock) and is written independently.
* nwlsm is free software: you can redistribute it and/or modify it under the terms of the GNU General
* Public License, version 3, as published by the Free Software Foundation. It is distributed in the hope
* that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
* or FITNESS FOR A PARTICULAR PURPOSE. See the file LICENSE for the full text.

capture program drop nwlsm
program nwlsm, eclass
	version 14
	gettoken sub rest : 0, parse(" ,")
	if "`sub'" == "version" {
		nwlsm_version
		exit
	}
	di as error "nwlsm is under development: no estimation command yet (see {bf:nwlsm version})"
	exit 199
end
