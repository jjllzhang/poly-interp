.PHONY: release debug test bench-interp bench-flint plot-interp plot-flint clean

release:
	cmake --preset release
	cmake --build --preset release --target interp_bench bench_flint

debug:
	cmake --preset debug
	cmake --build --preset debug

test: release
	ctest --preset release

bench-interp: release
	mkdir -p data
	./build/release/interp_bench --prime=all --min_pow=10 --max_pow=20 --repeats=3 --csv=1 > data/bench.csv

bench-flint: release
	mkdir -p data
	./build/release/bench_flint all --csv data/bench_flint.csv

plot-interp: release
	cmake --build --preset release --target plot_interp

plot-flint: release
	cmake --build --preset release --target plot_flint

clean:
	rm -rf build
