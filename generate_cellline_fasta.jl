proj_info = (id = "cov2",)
using Pkg
Pkg.activate(@__DIR__)

using Revise, FastaIO

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")

includet(joinpath(misc_scripts_path, "fasta_reader.jl"))
includet(joinpath(misc_scripts_path, "saav_utils.jl"))

proteins_df = Fasta.read_ensembl(joinpath(data_path, "msfasta/ENSEMBL_A549_HEK293_HeLa_with_SAAVs.fasta"), has_annotations=false, has_SAAVs=true)

fastaout = FastaIO.FastaWriter(joinpath(data_path, "msfasta/ENSEMBL_A549_HEK293_HeLa_with_SAAVs_20200812.fasta"))
for r in eachrow(proteins_df)
    SAAVUtils.writeentry(fastaout, r.protein_ac, r.seq, !ismissing(r.SAAVs) ? parse.(Ref(SAAVUtils.SAAV), split(r.SAAVs, ';')) : nothing)
end
close(fastaout)
