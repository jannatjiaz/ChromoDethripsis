#!/usr/bin/julia

ENV["GKSwstype"]="100"

using Bio
using BioAlignments
using BioSequences
using BioSymbols
using GenomicFeatures
using NamedArrays
using StatsBase
using FastaIO
using SpecialFunctions
using FileIO
using GeneticVariation
using DataFrames
using HypothesisTests
using Distributions
using FreqTables
using DelimitedFiles
using CSV
using Plots

function haplotype_assign(sv_file::String, hap_bam::String, het_snp_basects::String,
                            phased_het_SNPs::String, output_file::String, chr::String)

    # First, read in SV data and generate haplotype calls associated with each SV
    sv = open(sv_file, "r")
    hap_bam_reader = open(BAM.Reader, hap_bam, index= hap_bam * ".bai")
    sv_hap_dict = Dict{Int,Dict{Int,Int}}()

    for curr_sv in eachline(sv)
        if curr_sv[1] == '#'
            continue
        end
        sv_data = split(curr_sv, "\t")
        sv_readnames = split(sv_data[9], ";")
        if sv_data[1] == chr
            sv_hap_dict = extract_sv_haplotype_from_reads(hap_bam_reader, chr, parse(Int, sv_data[2]),
                                        sv_readnames, sv_hap_dict)
        end
        if sv_data[4] == chr
            sv_hap_dict = extract_sv_haplotype_from_reads(hap_bam_reader, chr, parse(Int, sv_data[5]),
                                        sv_readnames, sv_hap_dict)
        end
    end


end


"""
    extract_sv_haplotype_from_reads(hap_bam_reader, chr, pos1, sv_readnames, sv_hap_dict)

Internal function to extract haplotype information from reads in `sv_readnames` crossing `chr`:`pos1`
    from the WhatsHap phasing information contained in `hap_bam_reader`
"""
function extract_sv_haplotype_from_reads(hap_bam_reader, chr, pos1, sv_readnames, sv_hap_dict)
    for record in eachoverlap(hap_bam_reader, Interval(chr, pos1-1000, pos1+1000))
        if BAM.tempname(record) ∈ sv_readnames # Read spans the current SV
            if haskey(record, "PS")
                if haskey(sv_hap_dict, record["PS"])
                    if haskey(sv_hap_dict[record["PS"]], record["HP"])
                        sv_hap_dict[record["PS"]][record["HP"]] += 1
                    else
                        sv_hap_dict[record["PS"]][record["HP"]] = 1
                    end
                else
                    sv_hap_dict[record["PS"]] = Dict(record["HP"] => 1)
                end
            end
        end
    end
    return(sv_hap_dict)
end


"""
    phasing_het_SNPs(vcf_file::String, region::Interval, basects_file::String, plot_file::String)

Function to read in SNPs from `basects_file`, generated as output of `mpileup_parse.jl`, take phasing
    information from `vcf_file` and return (1) a `plot_file` of depth and VAF across genomic interval in
    `region` and (2) a DataFrame of basecounts for the genomic region with phasing information

The phasing information in the output is carried in 2 extra columns of the output:
    `Phase_set` => the name of the haplotype block the SNP falls in;
    `Phased_allele` =>  "0|1" means reference allele on haplotype block 1, while
                        "1|0" means reference allele on haplotype block 2.
"""
function phasing_het_SNPs(vcf_file::String, region::Interval, basects_file::String, plot_file::String)
    vcf_in = open(vcf_file, "r")
    basects = CSV.read(basects_file, DataFrame, delim = "\t")
    basects = basects[(basects.Chr .== seqname(region)) .&
                (leftposition(region) .<= basects.Pos .<= rightposition(region)), :]
    basects.Phase_set = repeat(["Unphased";], size(basects,1))
    basects.Phased_allele = repeat(["Unphased";], size(basects,1))
    phase_dict = Dict{String,Array{String,1}}()

    for line in eachline(vcf_in)
        if line[1] == '#' # Skip header
            continue
        end
        line_spl = split(line, "\t")
        pos = parse(Int, line_spl[2])
        if !(line_spl[1] == seqname(region) && leftposition(region) <= pos <= rightposition(region))
            continue # Outside genomic interval of interest
        end

        format_split = split(line_spl[9], ":")
        if "PS" ∉ format_split # SNP not phased to haplotype block
            continue
        end

        sample_split = split(line_spl[10], ":")
        format_dict = Dict(string.(format_split) .=> string.(sample_split))

        if format_dict["GT"] ∈ ["0|1","1|0"]
            phase_dict["$(line_spl[1]):$pos"] = [format_dict["PS"], format_dict["GT"]]
        end
    end

    for (ix, curr_row) in enumerate(eachrow(basects))
        if haskey(phase_dict, "$(curr_row.Chr):$(curr_row.Pos)")
            basects[ix,"Phase_set"] = phase_dict["$(curr_row.Chr):$(curr_row.Pos)"][1]
            basects[ix,"Phased_allele"] = phase_dict["$(curr_row.Chr):$(curr_row.Pos)"][2]
        end
    end

    basects.Depth = basects.Ref_count .+ basects.Alt_count
    basects.VAF = basects.Alt_count ./ basects.Depth
    pl = plot_depth_VAF(basects)
    savefig(pl, plot_file)
    return(basects)
end


"""
    function plot_depth_VAF(basects)

Function to generate plot of depth and VAF for SNP data included in `basects`
"""
function plot_depth_VAF(basects)
    pl1 = plot(basects.Pos[basects.Phased_allele .== "0|1"], basects.Depth[basects.Phased_allele .== "0|1"],
                seriestype = :scatter, seriesalpha = 0.5, markerstrokewidth=-1, markersize=2,
                markerstrokecolor=:blue2, markercolor=:blue2, xlabel="Genomic position", ylabel="Depth")
    plot!(basects.Pos[basects.Phased_allele .== "1|0"], basects.Depth[basects.Phased_allele .== "1|0"],
                seriestype = :scatter, seriesalpha = 0.5, markerstrokewidth=-1, markersize=2,
                markerstrokecolor=:red2, markercolor=:red2, legend=false)
    pl2 = plot(basects.Pos[basects.Phased_allele .== "0|1"], basects.VAF[basects.Phased_allele .== "0|1"],
                seriestype = :scatter, seriesalpha = 0.5, markerstrokewidth=-1, markersize=2,
                markerstrokecolor=:blue2, markercolor=:blue2, xlabel="Genomic position", ylabel="VAF")
    plot!(basects.Pos[basects.Phased_allele .== "1|0"], basects.VAF[basects.Phased_allele .== "1|0"],
                seriestype = :scatter, seriesalpha = 0.5, markerstrokewidth=-1, markersize=2,
                markerstrokecolor=:red2, markercolor=:red2, legend=false)
    return(plot(pl1, pl2, layout=(2,1)))
end


"""
    allele_fraction_hap_blocks(basects)

Function to take `basects` DataFrame from `phasing_het_SNPs()` and summarise the average read depths
    for alleles on block 1 and block 2 of each haplotype and average VAF
"""
function allele_fraction_hap_blocks(basects::DataFrame)
    hap_ct = freqtable(basects.Phase_set)
    output = DataFrame(Haplotype = names(hap_ct)[1], Number_SNPs = hap_ct)
    output.ave_VAF_block1 = zeros(size(output,1))
    output.ave_depth_block1 = zeros(size(output,1))
    output.ave_depth_block2 = zeros(size(output,1))

    for (ix,block) in enumerate(eachrow(output))
        if block.Haplotype != "Unphased"
            block1_bool = ((basects.Phase_set .== block.Haplotype) .& (basects.Phased_allele .== "0|1"))
            block2_bool = ((basects.Phase_set .== block.Haplotype) .& (basects.Phased_allele .== "1|0"))
            reads_block1 = sum(basects.Ref_count[block1_bool]) + sum(basects.Alt_count[block2_bool])
            reads_block2 = sum(basects.Alt_count[block1_bool]) + sum(basects.Ref_count[block2_bool])
            output.ave_depth_block1[ix] = reads_block1 / block.Number_SNPs
            output.ave_depth_block2[ix] = reads_block2 / block.Number_SNPs
            output.ave_VAF_block1[ix] = reads_block1 / (reads_block1 + reads_block2)
        end
    end

    output.Haplotype_int = [parse.(Int64, output.Haplotype[1:(end-1)]); 5000000000]
    sort!(output, :Haplotype_int)
    return(output[:,Not(:Haplotype_int)])
end

