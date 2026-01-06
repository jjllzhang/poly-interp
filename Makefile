.PHONY: release debug test bench-interp bench-ops bench-ops-flint bench-flint plot-interp plot-interp-compare plot-ops plot-ops-flint plot-ops-compare plot-flint clean

release:
	cmake --preset release
	cmake --build --preset release --target interp_bench ops_bench bench_flint ops_bench_flint

debug:
	cmake --preset debug
	cmake --build --preset debug

test: release
	ctest --preset release

bench-interp: release
	mkdir -p data
	./build/release/interp_bench --prime=all --min_pow=10 --max_pow=20 --repeats=3 --csv=1 > data/bench.csv

bench-ops: release
	mkdir -p data
	./build/release/ops_bench --prime=all --min_pow=8 --max_pow=20 --repeats=3 --csv=1 > data/bench_ops.csv

bench-ops-flint: release
	mkdir -p data
	./build/release/ops_bench_flint --prime=all --min_pow=8 --max_pow=20 --repeats=3 --csv=1 > data/bench_ops_flint.csv

plot-ops:
	mkdir -p plots/ops
	python3 scripts/plot.py data/bench_ops.csv plots/ops --format=ops --logy

plot-ops-flint:
	mkdir -p plots/ops_flint
	python3 scripts/plot.py data/bench_ops_flint.csv plots/ops_flint --format=ops --logy

plot-ops-compare:
	mkdir -p plots/ops_compare
	python3 scripts/plot.py data/bench_ops.csv plots/ops_compare --format=ops --compare-ops data/bench_ops_flint.csv --logy

bench-flint: release
	mkdir -p data
	./build/release/bench_flint all --csv data/bench_flint.csv

plot-interp: release
	cmake --build --preset release --target plot_interp

plot-interp-compare:
	mkdir -p plots/interp_compare
	python3 scripts/plot.py data/bench.csv plots/interp_compare --format=interp --compare-interp data/bench_flint.csv --logy

plot-flint: release
	cmake --build --preset release --target plot_flint

clean:
	rm -rf build
