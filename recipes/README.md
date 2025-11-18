
# Recipes for building RFDiffussionAA and LigandMPNN containers with Rocm support

Many (many) thanks to **Piero Coronica** of the Max Planck Computing and Data Facility who developed and provided these recipes.


Baisc tests using the examples provided with the githup repos are shown. Benchmarks on my system, using a single MI300A show good (>A100) performance.

## RFDiffusion

### BUILD IMAGE
```
apptainer build rf_diffusion_AA.sif rf_diffusion_AA.def
```

### RUN TEST

```
apptainer exec rf_diffusion_AA.sif python /rf_diffusion_all_atom/run_inference.py inference.deterministic=True diffuser.T=100 inference.output_prefix=output/ligand_only/sample inference.input_pdb=input/7v11.pdb contigmap.contigs=[\'150-150\'] inference.ligand=OQO inference.num_designs=1 inference.design_startnum=0
```

## LigandMPNN

### BUILD

```
apptainer build ligandmpnn.sif ligandmpnn.def
```

### RUN (from github README)

```
apptainer exec ligandmpnn.sif bash /LigandMPNN/get_model_params.sh "./model_params"
apptainer exec ligandmpnn.sif python /LigandMPNN/run.py --seed 111 --pdb_path "/LigandMPNN/inputs/1BC8.pdb" --out_folder "./outputs/default"
```
