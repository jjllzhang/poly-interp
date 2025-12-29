.PHONY: release debug test bench plot clean

release:
	cmake --preset release
	cmake --build --preset release

debug:
	cmake --preset debug
	cmake --build --preset debug

test: release
	ctest --preset release

bench: release
	cmake --build --preset release --target run_bench

plot: release
	cmake --build --preset release --target plot

clean:
	rm -rf build
