# -*- coding: utf-8 -*-
"""
rate_partitions.py
Written by: T Malm and C Peña - 20180517
Script to read a TIGER rate file (or other rate files) and produce partitions
from the rates usable for analyses.

Currently partition sizes are calculated as follows:
The first partition is calculated as

highest_rate - ((highest_rate - minimum_rate)/(divfactor)), the remaining bins as:
higher_rate - ((higher_rate - lower_rate)/divfactor + partition_number * 0.3)

leading to smaller bin rate spans with faster rates (=lower rates), the last
partition will be created when 10% or less characters are left.

Run the script from terminal with a rate-file, as first argument, and the wanted
divfactor, as second argument, and it will produce a new output file with
partition rate summaries and partition schemes for MrBayes and PHYLIP.

The rate file should be a file with rates separated by hard return (TIGER rate
output -rl command) and in order of the character alignment you want to use later.

Run as: "python rate_partitions.py X Y" where X is the input rate file and Y
is the division factor (e.g. 1.35)

The output file will be named after the input values: inputfilename_divfactor.txt
"""
import argparse
import os
import re
import sys


def run(infile, divfactor):
    """Executes computations and returns contents to be written in output file"""
    input_data = read_input_file(infile)
    infile_basename = os.path.basename(infile)
    
    max_rate, min_rate, num_chars = max(input_data), min(input_data), len(input_data)
    spread = max_rate - min_rate
    print("\nTotal data - {} sites.\n".format(num_chars))
    print("Slowest rate: {}".format(max_rate))
    print("Fastest rate: {}".format(min_rate))
    print("Rate spread: {}\n".format(spread))

    # setting values for partitioning
    upper_value = max_rate
    partitioned_sites_count = 0  # total partitioned sites count
    cutoff = num_chars * 0.1  # cutoff value for creating last partition
    output_phy = "PHYLIP  style"
    output_mrb = "MrBayes style\nbegin mrbayes;"
    output_info = make_output_description(
        divfactor, infile_basename, cutoff, max_rate, min_rate, num_chars, spread)

    bin_count = 0
    for partition in range(1, 100):
        number_remaining_sites = num_chars - partitioned_sites_count
        if number_remaining_sites <= cutoff:
            # for last partition to include all the rest
            lower_value = min_rate
            sites = generate_sites_last_partition(input_data, lower_value, upper_value)
        else:
            sites, lower_value = generate_sites(input_data, partition, upper_value, min_rate, divfactor)

        # info for output in file and screen
        output_info += "\nPartition_{}({} sites):	Rate-span: {}-{}\n".format(
            partition, len(sites), round(upper_value, 6), round(lower_value, 6))
        print("Partition_{} ({} sites): ".format(partition, len(sites)))
        print("Rate-span: {0:.5f}-{0:.5f}\n".format(upper_value, lower_value))

        output_mrb, output_phy = add_partitions_output(output_mrb, output_phy,
                                                       partition, sites)
        upper_value = lower_value  # resetting the upper range value to the current lower value (for next bin)
        partitioned_sites_count += len(sites)  # for total site count
        bin_count += 1

        # breaking loop on last partition
        if number_remaining_sites <= cutoff:
            break

    # more info
    if partitioned_sites_count != num_chars:
        msg = "Total sites partitioned is not identical to imported " \
              "sites!:{} vs {}".format(partitioned_sites_count, num_chars)
        print(msg)
        output_info += msg

    # fixing partition finishing for output for MrBayes partitioning
    partition_list = generate_partition_list(bin_count)
    out_finish = "\npartition Partitions = {}: {};".format(bin_count, partition_list)
    output_mrb += out_finish
    output_mrb += "\nset partition = Partitions;"

    # collecting outputs
    output_finished = [output_info, output_mrb, output_phy]
    output_finished = '\n\n\n'.join(output_finished)
    return output_finished


def generate_partition_list(bin_count):
    partitions = ""
    for partition in range(1, bin_count + 1):
        partition_string = "Partition_{}, ".format(partition)
        partitions += partition_string
    partitions = re.sub(", $", "", partitions)
    return partitions


def generate_sites(input_data, partition, upper_value, min_rate, divfactor):
    """Generate a list of sites whose evolutionary rate is between upper and lower values

    This is the most important function. This function contains all the logic
    to divide the characters in bins based on upper rate, minimum rate, partition
    number and divfactor.
    """
    if partition == 1:
        # for first partition
        lower_value = upper_value - ((upper_value - min_rate) / divfactor)
    else:
        # for all other partitions. other than the last one
        lower_value = upper_value - ((upper_value - min_rate) / (divfactor + partition * 0.3))

    sites = []
    for idx, rate in enumerate(input_data):
        if upper_value >= rate > lower_value:
            sites.append(idx + 1)
    return sites, lower_value


def generate_sites_last_partition(input_data, lower_value, upper_value):
    """Puts all remaining sites in a last bin"""
    sites = []
    for idx, rate in enumerate(input_data):
        if upper_value > rate >= lower_value:
            sites.append(idx + 1)
    return sites


def add_partitions_output(output_mrb, output_phy, partition, sites):
    """Adds to output the Partition line with partition number and list of characters"""
    if sites:
        # setting the output format for charsets as MrBayes partitions
        charset = " ".join([str(site) for site in sites])
        output_mrb += "\nCharset Partition_{} = {};".format(partition, charset)

        # setting the output for phylip partitions
        charset = ", ".join([str(site) for site in sites])
        output_phy += "\nDNA, Partition_{} = {}".format(partition, charset)

    return output_mrb, output_phy


def make_output_description(divfactor, infile_basename, cutoff, max_rate,
                            min_rate, num_chars, spread):
    """Generate text for output for section that describes method"""
    output = "Partition output from ratepartitions.py\n--Written by Tobias Malm " \
         "(20121130)\n\nFor rate file: {} with {} sites!\n".format(
        infile_basename, num_chars)
    output += "\nManually set dividing factor: {}\n".format(divfactor)
    output += "Partitions calculated according to:\n\t1st partition: highest " \
          "rate - ((highest rate - minimum_rate)/({0})),\n\tthe remaining as " \
          "lower boundary rate= upper boundary rate -((upper boundary rate - " \
          "minimum rate)/({0}+Partitionnumber*0.3)).\n".format(divfactor)
    output += "\tLast partition created when less than 10% of total characters " \
          "are left (={} characters).\n".format(cutoff)
    output += "\nRate spread of entire data set (Highest (slowest, 1=invariant) to " \
          "lowest (fastest) ): Highest: {}, lowest: {}, spread: {}".format(
        max_rate, min_rate, spread)
    return output


def read_input_file(infile):
    with open(infile, "r") as handle:
        lines = handle.readlines()

        # reading rates lines into List
        rate_values = []

        for line in lines:
            clean_line = line.strip()
            if clean_line:
                rate_values.append(float(clean_line))

        return rate_values


def write_output_file(output_data, infile, divfactor):
    # opening outfile and writing to it
    outfile = "{}_{}.txt".format(infile, divfactor)
    print("Output file with rate partition summary and MrBayes and PHYLIP "
          "partition schemes has been created: {}".format(outfile))

    with open(outfile, "w") as handle:
        handle.write(output_data)


def verify_divfactor(divfactor):
    """Check that divfactor is greater or equal than value 1.1"""
    error_msg = "You need to enter factor for division as positive numerical " \
                "value (greater or equal than 1.1)"
    if divfactor < 1.1:
        print(error_msg)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "rate_file_txt",
        help="A file with rates separated by hard return (TIGER rate output "
             "-rl command)",
    )
    parser.add_argument(
        "divfactor",
        type=float,
        help="Factor for division as positive numerical value (greater or "
             "equal than 1.1)",
    )

    args = parser.parse_args()
    infile = args.rate_file_txt
    divfactor = args.divfactor

    verify_divfactor(divfactor)

    output_data = run(infile, divfactor)
    write_output_file(output_data, infile, divfactor)


if __name__ == "__main__":
    main()
