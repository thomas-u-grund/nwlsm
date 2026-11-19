# Native (C) Stata plugins of the nwlsm package: lsm_sim.
# Output: ../lib/plugins/<platform>/<plugin>_<platform>.plugin (the platform suffix matters for `net install`,
# which flattens all files into one folder). spi/ holds StataCorp's plugin interface (SPI 3.0, redistributed as
# published at https://www.stata.com/plugins/); it is not this package's own code.
#   make macos    universal (x86_64 + arm64) binaries, run on macOS
#   make unix     Linux (run on a Linux host or the CI runner)
# Windows (Visual Studio x64 prompt): cl /LD /O2 /DSYSTEM=STWIN32 /I spi <plugin>.c spi\stplugin.c /Fe:..\lib\plugins\windows\<plugin>_windows.plugin

CC = clang
CFLAGS = -O2 -Wall -Ispi

.PHONY: macos unix clean macos-lsm_sim unix-lsm_sim

macos: macos-lsm_sim

unix: unix-lsm_sim

macos-lsm_sim:
	mkdir -p ../lib/plugins/macos
	$(CC) -bundle -o lsm_sim.x86_64 lsm_sim.c spi/stplugin.c \
		-DSYSTEM=APPLEMAC -target x86_64-apple-macos10.11 $(CFLAGS)
	$(CC) -bundle -o lsm_sim.arm64 lsm_sim.c spi/stplugin.c \
		-DSYSTEM=APPLEMAC -target arm64-apple-macos11 $(CFLAGS)
	lipo -create -output ../lib/plugins/macos/lsm_sim_macos.plugin lsm_sim.x86_64 lsm_sim.arm64
	codesign -s - ../lib/plugins/macos/lsm_sim_macos.plugin
	rm -f lsm_sim.x86_64 lsm_sim.arm64
	@echo "Built ../lib/plugins/macos/lsm_sim_macos.plugin"

unix-lsm_sim:
	mkdir -p ../lib/plugins/unix
	gcc -shared -fPIC -DSYSTEM=OPUNIX $(CFLAGS) -pthread -o ../lib/plugins/unix/lsm_sim_unix.plugin \
		lsm_sim.c spi/stplugin.c -ldl -lpthread
	@echo "Built ../lib/plugins/unix/lsm_sim_unix.plugin"

clean:
	rm -f lsm_sim.x86_64 lsm_sim.arm64 ../lib/plugins/macos/lsm_sim_macos.plugin ../lib/plugins/unix/lsm_sim_unix.plugin
