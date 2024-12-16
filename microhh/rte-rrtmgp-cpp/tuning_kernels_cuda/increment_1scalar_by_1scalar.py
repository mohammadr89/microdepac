import argparse
import json
from collections import OrderedDict
import kernel_tuner as kt
import numpy as np

import common


# Parse command line arguments
def parse_command_line():
    parser = argparse.ArgumentParser(description="Tuning script for increment_1scalar_by_1scalar kernels")
    parser.add_argument("--tune", default=False, action="store_true")
    parser.add_argument("--run", default=False, action="store_true")
    parser.add_argument("--best_configuration", default=False, action="store_true")
    parser.add_argument("--ncol", type=int, default=512)
    parser.add_argument("--nlay", type=int, default=140)
    parser.add_argument("--ngpt", type=int, default=224)
    parser.add_argument("--nbnd", type=int, default=14)
    parser.add_argument("--block_size_x", type=int, default=32)
    parser.add_argument("--block_size_y", type=int, default=16)
    parser.add_argument("--block_size_z", type=int, default=1)
    parser.add_argument("--block_size_x_bybnd", type=int, default=16)
    parser.add_argument("--block_size_y_bybnd", type=int, default=4)
    parser.add_argument("--block_size_z_bybnd", type=int, default=1)
    parser.add_argument("--loop_unroll_bybnd", type=int, default=1)
    return parser.parse_args()


# Run one instance of the kernel and test output
def run_and_test(params: OrderedDict):
    global tau1
    global tau1_ref
    # kernel
    print(f"Running {kernel_name['kernel']} [{params['kernel']['block_size_x']}, {params['kernel']['block_size_y']}, {params['kernel']['block_size_z']}]")
    increment_1scalar_by_1scalar(tau1_ref)
    result = kt.run_kernel(kernel_name["kernel"], kernels_src, problem_size, args["kernel"], params["kernel"],
                           compiler_options=common.cp)
    common.compare_fields(tau1_ref, result[3], "tau1")
    # bybnd
    print(
        f"Running {kernel_name['bybnd']} [{params['bybnd']['block_size_x']}, {params['bybnd']['block_size_y']}, {params['bybnd']['block_size_z']}, {params['bybnd']['loop_unroll_factor_nbnd']}]")
    tau1 = common.random(size, common.type_float)
    tau1_ref = tau1
    args["bybnd"][3] = tau1
    args["bybnd_ref"][3] = tau1_ref
    ref_result = kt.run_kernel(kernel_name["bybnd"], kernels_src_ref, problem_size, args["bybnd_ref"],
                               OrderedDict(block_size_x=16, block_size_y=4, block_size_z=1), compiler_options=common.cp)
    result = kt.run_kernel(kernel_name["bybnd"], kernels_src, problem_size, args["bybnd"], params["bybnd"],
                           compiler_options=common.cp)
    common.compare_fields(ref_result[3], result[3], "tau1")


# Tuning
def tune():
    global tau1
    global tau1_ref
    tune_params = OrderedDict()
    tune_params["block_size_x"] = [2 ** i for i in range(0, 11)]
    tune_params["block_size_y"] = [2 ** i for i in range(0, 11)]
    tune_params["block_size_z"] = [2 ** i for i in range(0, 7)]
    restrictions = [f"block_size_x <= {ncol}", f"block_size_y <= {nlay}", f"block_size_z <= {ngpt}"]
    # kernel
    print(f"Tuning {kernel_name['kernel']}")
    increment_1scalar_by_1scalar(tau1_ref)
    answer = [None for _ in range(0, len(args["kernel"]))]
    answer[3] = tau1_ref
    result, env = kt.tune_kernel(kernel_name["kernel"], kernels_src, problem_size, args["kernel"], tune_params,
                                 answer=answer, compiler_options=common.cp, verbose=True, restrictions=restrictions)
    with open("timings_increment_1scalar_by_1scalar_kernel.json", "w") as fp:
        json.dump(result, fp)
    # bybnd
    print(f"Tuning {kernel_name['bybnd']}")
    tau1 = common.random(size, common.type_float)
    tau1_ref = tau1
    args["bybnd"][3] = tau1
    args["bybnd_ref"][3] = tau1_ref
    ref_result = kt.run_kernel(kernel_name["bybnd"], kernels_src_ref, problem_size, args["bybnd_ref"],
                               OrderedDict(block_size_x=16, block_size_y=4, block_size_z=1), compiler_options=common.cp)
    answer = [None for _ in range(0, len(args["bybnd"]))]
    answer[3] = ref_result[3]
    tune_params["loop_unroll_factor_nbnd"] = [0] + [i for i in range(1, nbnd + 1) if nbnd // i == nbnd / i]
    result, env = kt.tune_kernel(kernel_name["bybnd"], kernels_src, problem_size, args["bybnd"], tune_params,
                                 answer=answer, compiler_options=common.cp, verbose=True, restrictions=restrictions)
    with open("timings_inc_1scalar_by_1scalar_bybnd_kernel.json", "w") as fp:
        json.dump(result, fp)


# Python implementation, for correctness only
def increment_1scalar_by_1scalar(out: np.ndarray):
    for icol in range(0, ncol):
        for ilay in range(0, nlay):
            for igpt in range(0, ngpt):
                index = (igpt * nlay * ncol) + (ilay * ncol) + icol
                out[index] = tau2[index] * 2.0


if __name__ == "__main__":
    command_line = parse_command_line()

    kernels_src = common.dir_name + "../src_kernels_cuda/optical_props_kernels.cu"
    kernels_src_ref = common.dir_name + "/reference_kernels/optical_props_kernels.cu"

    # Data
    ncol = common.type_int(command_line.ncol)
    nlay = common.type_int(command_line.nlay)
    ngpt = common.type_int(command_line.ngpt)
    nbnd = common.type_int(command_line.nbnd)
    size = ncol * nlay * ngpt
    tau1 = common.random(size, common.type_float)
    tau1_ref = tau1
    tau2 = common.random(size, common.type_float)
    band_lims_gpt = np.fromfile(f"{common.bin_path}/band_lims_gpt.bin", dtype=common.type_int)

    kernel_name = OrderedDict()
    kernel_name["kernel"] = f"increment_1scalar_by_1scalar_kernel<{common.str_float}>"
    kernel_name["bybnd"] = f"inc_1scalar_by_1scalar_bybnd_kernel<{common.str_float}>"
    problem_size = (ncol, nlay, ngpt)
    args = OrderedDict()
    args["kernel"] = [ncol, nlay, ngpt, tau1, tau2]
    args["bybnd"] = [ncol, nlay, ngpt, tau1, tau2, nbnd, band_lims_gpt]
    args["bybnd_ref"] = [ncol, nlay, ngpt, tau1_ref, tau2, nbnd, band_lims_gpt]

    if command_line.tune:
        tune()
    elif command_line.run:
        parameters = OrderedDict()
        parameters["kernel"] = OrderedDict()
        parameters["bybnd"] = OrderedDict()
        if command_line.best_configuration:
            best_configuration = common.best_configuration("timings_increment_1scalar_by_1scalar_kernel.json")
            parameters["kernel"]["block_size_x"] = best_configuration["block_size_x"]
            parameters["kernel"]["block_size_y"] = best_configuration["block_size_y"]
            parameters["kernel"]["block_size_z"] = best_configuration["block_size_z"]
            best_configuration = common.best_configuration("timings_inc_1scalar_by_1scalar_bybnd_kernel.json")
            parameters["bybnd"]["block_size_x"] = best_configuration["block_size_x"]
            parameters["bybnd"]["block_size_y"] = best_configuration["block_size_y"]
            parameters["bybnd"]["block_size_z"] = best_configuration["block_size_z"]
            parameters["bybnd"]["loop_unroll_factor_nbnd"] = best_configuration["loop_unroll_factor_nbnd"]
        else:
            parameters["kernel"]["block_size_x"] = command_line.block_size_x
            parameters["kernel"]["block_size_y"] = command_line.block_size_y
            parameters["kernel"]["block_size_z"] = command_line.block_size_z
            parameters["bybnd"]["block_size_x"] = command_line.block_size_x_bybnd
            parameters["bybnd"]["block_size_y"] = command_line.block_size_y_bybnd
            parameters["bybnd"]["block_size_z"] = command_line.block_size_z_bybnd
            parameters["bybnd"]["loop_unroll_factor_nbnd"] = command_line.loop_unroll_bybnd
        run_and_test(parameters)
